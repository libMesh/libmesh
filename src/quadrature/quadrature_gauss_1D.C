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



// C++ includes

// Local includes
#include "libmesh/quadrature_gauss.h"

namespace libMesh
{



void QGauss::init_1D(const ElemType, unsigned int)
{
  //----------------------------------------------------------------------
  // 1D quadrature rules
  switch(get_order())
    {
    case CONSTANT:
    case FIRST:
      {
        _points.resize (1);
        _weights.resize(1);

        _points[0](0)  = 0.;

        _weights[0]    = 2.;

        return;
      }
    case SECOND:
    case THIRD:
      {
        _points.resize (2);
        _weights.resize(2);

        _points[0](0) = Real(-5.7735026918962576450914878050196e-01L); // -sqrt(3)/3
        _points[1]    = -_points[0];

        _weights[0]   = 1.;
        _weights[1]   = _weights[0];

        return;
      }
    case FOURTH:
    case FIFTH:
      {
        _points.resize (3);
        _weights.resize(3);

        _points[ 0](0) = Real(-7.7459666924148337703585307995648e-01L);
        _points[ 1](0) = 0.;
        _points[ 2]    = -_points[0];

        _weights[ 0]   = Real(5.5555555555555555555555555555556e-01L);
        _weights[ 1]   = Real(8.8888888888888888888888888888889e-01L);
        _weights[ 2]   = _weights[0];

        return;
      }
    case SIXTH:
    case SEVENTH:
      {
        _points.resize (4);
        _weights.resize(4);

        _points[ 0](0) = Real(-8.6113631159405257522394648889281e-01L);
        _points[ 1](0) = Real(-3.3998104358485626480266575910324e-01L);
        _points[ 2]    = -_points[1];
        _points[ 3]    = -_points[0];

        _weights[ 0]   = Real(3.4785484513745385737306394922200e-01L);
        _weights[ 1]   = Real(6.5214515486254614262693605077800e-01L);
        _weights[ 2]   = _weights[1];
        _weights[ 3]   = _weights[0];

        return;
      }
    case EIGHTH:
    case NINTH:
      {
        _points.resize (5);
        _weights.resize(5);

        _points[ 0](0) = Real(-9.0617984593866399279762687829939e-01L);
        _points[ 1](0) = Real(-5.3846931010568309103631442070021e-01L);
        _points[ 2](0) = 0.;
        _points[ 3]    = -_points[1];
        _points[ 4]    = -_points[0];

        _weights[ 0]   = Real(2.3692688505618908751426404071992e-01L);
        _weights[ 1]   = Real(4.7862867049936646804129151483564e-01L);
        _weights[ 2]   = Real(5.6888888888888888888888888888889e-01L);
        _weights[ 3]   = _weights[1];
        _weights[ 4]   = _weights[0];

        return;
      }
    case TENTH:
    case ELEVENTH:
      {
        _points.resize (6);
        _weights.resize(6);

        _points[ 0](0) = Real(-9.3246951420315202781230155449399e-01L);
        _points[ 1](0) = Real(-6.6120938646626451366139959501991e-01L);
        _points[ 2](0) = Real(-2.3861918608319690863050172168071e-01L);
        _points[ 3]    = -_points[2];
        _points[ 4]    = -_points[1];
        _points[ 5]    = -_points[0];

        _weights[ 0]   = Real(1.7132449237917034504029614217273e-01L);
        _weights[ 1]   = Real(3.6076157304813860756983351383772e-01L);
        _weights[ 2]   = Real(4.6791393457269104738987034398955e-01L);
        _weights[ 3]   = _weights[2];
        _weights[ 4]   = _weights[1];
        _weights[ 5]   = _weights[0];

        return;
      }
    case TWELFTH:
    case THIRTEENTH:
      {
        _points.resize (7);
        _weights.resize(7);

        _points[ 0](0) = Real(-9.4910791234275852452618968404785e-01L);
        _points[ 1](0) = Real(-7.4153118559939443986386477328079e-01L);
        _points[ 2](0) = Real(-4.0584515137739716690660641207696e-01L);
        _points[ 3](0) = 0.;
        _points[ 4]    = -_points[2];
        _points[ 5]    = -_points[1];
        _points[ 6]    = -_points[0];

        _weights[ 0]   = Real(1.2948496616886969327061143267908e-01L);
        _weights[ 1]   = Real(2.7970539148927666790146777142378e-01L);
        _weights[ 2]   = Real(3.8183005050511894495036977548898e-01L);
        _weights[ 3]   = Real(4.1795918367346938775510204081633e-01L);
        _weights[ 4]   = _weights[2];
        _weights[ 5]   = _weights[1];
        _weights[ 6]   = _weights[0];

        return;
      }
    case FOURTEENTH:
    case FIFTEENTH:
      {
        _points.resize (8);
        _weights.resize(8);

        _points[ 0](0) = Real(-9.6028985649753623168356086856947e-01L);
        _points[ 1](0) = Real(-7.9666647741362673959155393647583e-01L);
        _points[ 2](0) = Real(-5.2553240991632898581773904918925e-01L);
        _points[ 3](0) = Real(-1.8343464249564980493947614236018e-01L);
        _points[ 4]    = -_points[3];
        _points[ 5]    = -_points[2];
        _points[ 6]    = -_points[1];
        _points[ 7]    = -_points[0];

        _weights[ 0]   = Real(1.0122853629037625915253135430996e-01L);
        _weights[ 1]   = Real(2.2238103445337447054435599442624e-01L);
        _weights[ 2]   = Real(3.1370664587788728733796220198660e-01L);
        _weights[ 3]   = Real(3.6268378337836198296515044927720e-01L);
        _weights[ 4]   = _weights[3];
        _weights[ 5]   = _weights[2];
        _weights[ 6]   = _weights[1];
        _weights[ 7]   = _weights[0];

        return;
      }
    case SIXTEENTH:
    case SEVENTEENTH:
      {
        _points.resize (9);
        _weights.resize(9);

        _points[ 0](0) = Real(-9.6816023950762608983557620290367e-01L);
        _points[ 1](0) = Real(-8.3603110732663579429942978806973e-01L);
        _points[ 2](0) = Real(-6.1337143270059039730870203934147e-01L);
        _points[ 3](0) = Real(-3.2425342340380892903853801464334e-01L);
        _points[ 4](0) = 0.;
        _points[ 5]    = -_points[3];
        _points[ 6]    = -_points[2];
        _points[ 7]    = -_points[1];
        _points[ 8]    = -_points[0];

        _weights[ 0]   = Real(8.1274388361574411971892158110524e-02L);
        _weights[ 1]   = Real(1.8064816069485740405847203124291e-01L);
        _weights[ 2]   = Real(2.6061069640293546231874286941863e-01L);
        _weights[ 3]   = Real(3.1234707704000284006863040658444e-01L);
        _weights[ 4]   = Real(3.3023935500125976316452506928697e-01L);
        _weights[ 5]   = _weights[3];
        _weights[ 6]   = _weights[2];
        _weights[ 7]   = _weights[1];
        _weights[ 8]   = _weights[0];

        return;
      }
    case EIGHTTEENTH:
    case NINETEENTH:
      {
        _points.resize (10);
        _weights.resize(10);

        _points[ 0](0) = Real(-9.7390652851717172007796401208445e-01L);
        _points[ 1](0) = Real(-8.6506336668898451073209668842349e-01L);
        _points[ 2](0) = Real(-6.7940956829902440623432736511487e-01L);
        _points[ 3](0) = Real(-4.3339539412924719079926594316578e-01L);
        _points[ 4](0) = Real(-1.4887433898163121088482600112972e-01L);
        _points[ 5]    = -_points[4];
        _points[ 6]    = -_points[3];
        _points[ 7]    = -_points[2];
        _points[ 8]    = -_points[1];
        _points[ 9]    = -_points[0];

        _weights[ 0]   = Real(6.6671344308688137593568809893332e-02L);
        _weights[ 1]   = Real(1.4945134915058059314577633965770e-01L);
        _weights[ 2]   = Real(2.1908636251598204399553493422816e-01L);
        _weights[ 3]   = Real(2.6926671930999635509122692156947e-01L);
        _weights[ 4]   = Real(2.9552422471475287017389299465134e-01L);
        _weights[ 5]   = _weights[4];
        _weights[ 6]   = _weights[3];
        _weights[ 7]   = _weights[2];
        _weights[ 8]   = _weights[1];
        _weights[ 9]   = _weights[0];

        return;
      }

    case TWENTIETH:
    case TWENTYFIRST:
      {
        _points.resize (11);
        _weights.resize(11);

        _points[ 0](0) = Real(-9.7822865814605699280393800112286e-01L);
        _points[ 1](0) = Real(-8.8706259976809529907515776930393e-01L);
        _points[ 2](0) = Real(-7.3015200557404932409341625203115e-01L);
        _points[ 3](0) = Real(-5.1909612920681181592572566945861e-01L);
        _points[ 4](0) = Real(-2.6954315595234497233153198540086e-01L);
        _points[ 5](0) = 0.;
        _points[ 6]    = -_points[4];
        _points[ 7]    = -_points[3];
        _points[ 8]    = -_points[2];
        _points[ 9]    = -_points[1];
        _points[10]    = -_points[0];

        _weights[ 0]   = Real(5.5668567116173666482753720442549e-02L);
        _weights[ 1]   = Real(1.2558036946490462463469429922394e-01L);
        _weights[ 2]   = Real(1.8629021092773425142609764143166e-01L);
        _weights[ 3]   = Real(2.3319376459199047991852370484318e-01L);
        _weights[ 4]   = Real(2.6280454451024666218068886989051e-01L);
        _weights[ 5]   = Real(2.7292508677790063071448352833634e-01L);
        _weights[ 6]   = _weights[4];
        _weights[ 7]   = _weights[3];
        _weights[ 8]   = _weights[2];
        _weights[ 9]   = _weights[1];
        _weights[10]   = _weights[0];

        return;
      }

    case TWENTYSECOND:
    case TWENTYTHIRD:
      {
        _points.resize (12);
        _weights.resize(12);

        _points[ 0](0) = Real(-9.8156063424671925069054909014928e-01L);
        _points[ 1](0) = Real(-9.0411725637047485667846586611910e-01L);
        _points[ 2](0) = Real(-7.6990267419430468703689383321282e-01L);
        _points[ 3](0) = Real(-5.8731795428661744729670241894053e-01L);
        _points[ 4](0) = Real(-3.6783149899818019375269153664372e-01L);
        _points[ 5](0) = Real(-1.2523340851146891547244136946385e-01L);
        _points[ 6]    = -_points[5];
        _points[ 7]    = -_points[4];
        _points[ 8]    = -_points[3];
        _points[ 9]    = -_points[2];
        _points[10]    = -_points[1];
        _points[11]    = -_points[0];

        _weights[ 0]   = Real(4.7175336386511827194615961485017e-02L);
        _weights[ 1]   = Real(1.0693932599531843096025471819400e-01L);
        _weights[ 2]   = Real(1.6007832854334622633465252954336e-01L);
        _weights[ 3]   = Real(2.0316742672306592174906445580980e-01L);
        _weights[ 4]   = Real(2.3349253653835480876084989892483e-01L);
        _weights[ 5]   = Real(2.4914704581340278500056243604295e-01L);
        _weights[ 6]   = _weights[5];
        _weights[ 7]   = _weights[4];
        _weights[ 8]   = _weights[3];
        _weights[ 9]   = _weights[2];
        _weights[10]   = _weights[1];
        _weights[11]   = _weights[0];

        return;
      }

    case TWENTYFOURTH:
    case TWENTYFIFTH:
      {
        _points.resize (13);
        _weights.resize(13);

        _points[ 0](0) = Real(-9.8418305471858814947282944880711e-01L);
        _points[ 1](0) = Real(-9.1759839922297796520654783650072e-01L);
        _points[ 2](0) = Real(-8.0157809073330991279420648958286e-01L);
        _points[ 3](0) = Real(-6.4234933944034022064398460699552e-01L);
        _points[ 4](0) = Real(-4.4849275103644685287791285212764e-01L);
        _points[ 5](0) = Real(-2.3045831595513479406552812109799e-01L);
        _points[ 6](0) = 0.;
        _points[ 7]    = -_points[5];
        _points[ 8]    = -_points[4];
        _points[ 9]    = -_points[3];
        _points[10]    = -_points[2];
        _points[11]    = -_points[1];
        _points[12]    = -_points[0];

        _weights[ 0]   = Real(4.0484004765315879520021592200986e-02L);
        _weights[ 1]   = Real(9.2121499837728447914421775953797e-02L);
        _weights[ 2]   = Real(1.3887351021978723846360177686887e-01L);
        _weights[ 3]   = Real(1.7814598076194573828004669199610e-01L);
        _weights[ 4]   = Real(2.0781604753688850231252321930605e-01L);
        _weights[ 5]   = Real(2.2628318026289723841209018603978e-01L);
        _weights[ 6]   = Real(2.3255155323087391019458951526884e-01L);
        _weights[ 7]   = _weights[5];
        _weights[ 8]   = _weights[4];
        _weights[ 9]   = _weights[3];
        _weights[10]   = _weights[2];
        _weights[11]   = _weights[1];
        _weights[12]   = _weights[0];

        return;
      }

    case TWENTYSIXTH:
    case TWENTYSEVENTH:
      {
        _points.resize (14);
        _weights.resize(14);

        _points[ 0](0) = Real(-9.8628380869681233884159726670405e-01L);
        _points[ 1](0) = Real(-9.2843488366357351733639113937787e-01L);
        _points[ 2](0) = Real(-8.2720131506976499318979474265039e-01L);
        _points[ 3](0) = Real(-6.8729290481168547014801980301933e-01L);
        _points[ 4](0) = Real(-5.1524863635815409196529071855119e-01L);
        _points[ 5](0) = Real(-3.1911236892788976043567182416848e-01L);
        _points[ 6](0) = Real(-1.0805494870734366206624465021983e-01L);
        _points[ 7]    = -_points[6];
        _points[ 8]    = -_points[5];
        _points[ 9]    = -_points[4];
        _points[10]    = -_points[3];
        _points[11]    = -_points[2];
        _points[12]    = -_points[1];
        _points[13]    = -_points[0];

        _weights[ 0]   = Real(3.5119460331751863031832876138192e-02L);
        _weights[ 1]   = Real(8.0158087159760209805633277062854e-02L);
        _weights[ 2]   = Real(1.2151857068790318468941480907248e-01L);
        _weights[ 3]   = Real(1.5720316715819353456960193862384e-01L);
        _weights[ 4]   = Real(1.8553839747793781374171659012516e-01L);
        _weights[ 5]   = Real(2.0519846372129560396592406566122e-01L);
        _weights[ 6]   = Real(2.1526385346315779019587644331626e-01L);
        _weights[ 7]   = _weights[6];
        _weights[ 8]   = _weights[5];
        _weights[ 9]   = _weights[4];
        _weights[10]   = _weights[3];
        _weights[11]   = _weights[2];
        _weights[12]   = _weights[1];
        _weights[13]   = _weights[0];

        return;
      }

    case TWENTYEIGHTH:
    case TWENTYNINTH:
      {
        _points.resize (15);
        _weights.resize(15);

        _points[ 0](0) = Real(-9.8799251802048542848956571858661e-01L);
        _points[ 1](0) = Real(-9.3727339240070590430775894771021e-01L);
        _points[ 2](0) = Real(-8.4820658341042721620064832077422e-01L);
        _points[ 3](0) = Real(-7.2441773136017004741618605461394e-01L);
        _points[ 4](0) = Real(-5.7097217260853884753722673725391e-01L);
        _points[ 5](0) = Real(-3.9415134707756336989720737098105e-01L);
        _points[ 6](0) = Real(-2.0119409399743452230062830339460e-01L);
        _points[ 7](0) = 0.;
        _points[ 8]    = -_points[6];
        _points[ 9]    = -_points[5];
        _points[10]    = -_points[4];
        _points[11]    = -_points[3];
        _points[12]    = -_points[2];
        _points[13]    = -_points[1];
        _points[14]    = -_points[0];

        _weights[ 0]   = Real(3.0753241996117268354628393577204e-02L);
        _weights[ 1]   = Real(7.0366047488108124709267416450667e-02L);
        _weights[ 2]   = Real(1.0715922046717193501186954668587e-01L);
        _weights[ 3]   = Real(1.3957067792615431444780479451103e-01L);
        _weights[ 4]   = Real(1.6626920581699393355320086048121e-01L);
        _weights[ 5]   = Real(1.8616100001556221102680056186642e-01L);
        _weights[ 6]   = Real(1.9843148532711157645611832644384e-01L);
        _weights[ 7]   = Real(2.0257824192556127288062019996752e-01L);
        _weights[ 8]   = _weights[6];
        _weights[ 9]   = _weights[5];
        _weights[10]   = _weights[4];
        _weights[11]   = _weights[3];
        _weights[12]   = _weights[2];
        _weights[13]   = _weights[1];
        _weights[14]   = _weights[0];

        return;
      }

    case THIRTIETH:
    case THIRTYFIRST:
      {
        _points.resize (16);
        _weights.resize(16);

        _points[ 0](0) = Real(-9.8940093499164993259615417345033e-01L);
        _points[ 1](0) = Real(-9.4457502307323257607798841553461e-01L);
        _points[ 2](0) = Real(-8.6563120238783174388046789771239e-01L);
        _points[ 3](0) = Real(-7.5540440835500303389510119484744e-01L);
        _points[ 4](0) = Real(-6.1787624440264374844667176404879e-01L);
        _points[ 5](0) = Real(-4.5801677765722738634241944298358e-01L);
        _points[ 6](0) = Real(-2.8160355077925891323046050146050e-01L);
        _points[ 7](0) = Real(-9.5012509837637440185319335424958e-02L);
        _points[ 8]    = -_points[7];
        _points[ 9]    = -_points[6];
        _points[10]    = -_points[5];
        _points[11]    = -_points[4];
        _points[12]    = -_points[3];
        _points[13]    = -_points[2];
        _points[14]    = -_points[1];
        _points[15]    = -_points[0];

        _weights[ 0]   = Real(2.7152459411754094851780572456018e-02L);
        _weights[ 1]   = Real(6.2253523938647892862843836994378e-02L);
        _weights[ 2]   = Real(9.5158511682492784809925107602246e-02L);
        _weights[ 3]   = Real(1.2462897125553387205247628219202e-01L);
        _weights[ 4]   = Real(1.4959598881657673208150173054748e-01L);
        _weights[ 5]   = Real(1.6915651939500253818931207903033e-01L);
        _weights[ 6]   = Real(1.8260341504492358886676366796922e-01L);
        _weights[ 7]   = Real(1.8945061045506849628539672320828e-01L);
        _weights[ 8]   = _weights[7];
        _weights[ 9]   = _weights[6];
        _weights[10]   = _weights[5];
        _weights[11]   = _weights[4];
        _weights[12]   = _weights[3];
        _weights[13]   = _weights[2];
        _weights[14]   = _weights[1];
        _weights[15]   = _weights[0];

        return;
      }

    case THIRTYSECOND:
    case THIRTYTHIRD:
      {
        _points.resize (17);
        _weights.resize(17);

        _points[ 0](0) = Real(-9.9057547531441733567543401994067e-01L);
        _points[ 1](0) = Real(-9.5067552176876776122271695789580e-01L);
        _points[ 2](0) = Real(-8.8023915372698590212295569448816e-01L);
        _points[ 3](0) = Real(-7.8151400389680140692523005552048e-01L);
        _points[ 4](0) = Real(-6.5767115921669076585030221664300e-01L);
        _points[ 5](0) = Real(-5.1269053708647696788624656862955e-01L);
        _points[ 6](0) = Real(-3.5123176345387631529718551709535e-01L);
        _points[ 7](0) = Real(-1.7848418149584785585067749365407e-01L);
        _points[ 8](0) = 0.;
        _points[ 9]    = -_points[7];
        _points[10]    = -_points[6];
        _points[11]    = -_points[5];
        _points[12]    = -_points[4];
        _points[13]    = -_points[3];
        _points[14]    = -_points[2];
        _points[15]    = -_points[1];
        _points[16]    = -_points[0];

        _weights[ 0]   = Real(2.4148302868547931960110026287565e-02L);
        _weights[ 1]   = Real(5.5459529373987201129440165358245e-02L);
        _weights[ 2]   = Real(8.5036148317179180883535370191062e-02L);
        _weights[ 3]   = Real(1.1188384719340397109478838562636e-01L);
        _weights[ 4]   = Real(1.3513636846852547328631998170235e-01L);
        _weights[ 5]   = Real(1.5404576107681028808143159480196e-01L);
        _weights[ 6]   = Real(1.6800410215645004450997066378832e-01L);
        _weights[ 7]   = Real(1.7656270536699264632527099011320e-01L);
        _weights[ 8]   = Real(1.7944647035620652545826564426189e-01L);
        _weights[ 9]   = _weights[7];
        _weights[10]   = _weights[6];
        _weights[11]   = _weights[5];
        _weights[12]   = _weights[4];
        _weights[13]   = _weights[3];
        _weights[14]   = _weights[2];
        _weights[15]   = _weights[1];
        _weights[16]   = _weights[0];

        return;
      }

    case THIRTYFOURTH:
    case THIRTYFIFTH:
      {
        _points.resize (18);
        _weights.resize(18);

        _points[ 0](0) = Real(-9.9156516842093094673001600470615e-01L);
        _points[ 1](0) = Real(-9.5582394957139775518119589292978e-01L);
        _points[ 2](0) = Real(-8.9260246649755573920606059112715e-01L);
        _points[ 3](0) = Real(-8.0370495897252311568241745501459e-01L);
        _points[ 4](0) = Real(-6.9168704306035320787489108128885e-01L);
        _points[ 5](0) = Real(-5.5977083107394753460787154852533e-01L);
        _points[ 6](0) = Real(-4.1175116146284264603593179383305e-01L);
        _points[ 7](0) = Real(-2.5188622569150550958897285487791e-01L);
        _points[ 8](0) = Real(-8.4775013041735301242261852935784e-02L);
        _points[ 9]    = -_points[8];
        _points[10]    = -_points[7];
        _points[11]    = -_points[6];
        _points[12]    = -_points[5];
        _points[13]    = -_points[4];
        _points[14]    = -_points[3];
        _points[15]    = -_points[2];
        _points[16]    = -_points[1];
        _points[17]    = -_points[0];

        _weights[ 0]   = Real(2.1616013526483310313342710266452e-02L);
        _weights[ 1]   = Real(4.9714548894969796453334946202639e-02L);
        _weights[ 2]   = Real(7.6425730254889056529129677616637e-02L);
        _weights[ 3]   = Real(1.0094204410628716556281398492483e-01L);
        _weights[ 4]   = Real(1.2255520671147846018451912680020e-01L);
        _weights[ 5]   = Real(1.4064291467065065120473130375195e-01L);
        _weights[ 6]   = Real(1.5468467512626524492541800383637e-01L);
        _weights[ 7]   = Real(1.6427648374583272298605377646593e-01L);
        _weights[ 8]   = Real(1.6914238296314359184065647013499e-01L);
        _weights[ 9]   = _weights[8];
        _weights[10]   = _weights[7];
        _weights[11]   = _weights[6];
        _weights[12]   = _weights[5];
        _weights[13]   = _weights[4];
        _weights[14]   = _weights[3];
        _weights[15]   = _weights[2];
        _weights[16]   = _weights[1];
        _weights[17]   = _weights[0];

        return;
      }

    case THIRTYSIXTH:
    case THIRTYSEVENTH:
      {
        _points.resize (19);
        _weights.resize(19);

        _points[ 0](0) = Real(-9.9240684384358440318901767025326e-01L);
        _points[ 1](0) = Real(-9.6020815213483003085277884068765e-01L);
        _points[ 2](0) = Real(-9.0315590361481790164266092853231e-01L);
        _points[ 3](0) = Real(-8.2271465653714282497892248671271e-01L);
        _points[ 4](0) = Real(-7.2096617733522937861709586082378e-01L);
        _points[ 5](0) = Real(-6.0054530466168102346963816494624e-01L);
        _points[ 6](0) = Real(-4.6457074137596094571726714810410e-01L);
        _points[ 7](0) = Real(-3.1656409996362983199011732884984e-01L);
        _points[ 8](0) = Real(-1.6035864564022537586809611574074e-01L);
        _points[ 9](0) = 0.;
        _points[10]    = -_points[8];
        _points[11]    = -_points[7];
        _points[12]    = -_points[6];
        _points[13]    = -_points[5];
        _points[14]    = -_points[4];
        _points[15]    = -_points[3];
        _points[16]    = -_points[2];
        _points[17]    = -_points[1];
        _points[18]    = -_points[0];

        _weights[ 0]   = Real(1.9461788229726477036312041464438e-02L);
        _weights[ 1]   = Real(4.4814226765699600332838157401994e-02L);
        _weights[ 2]   = Real(6.9044542737641226580708258006013e-02L);
        _weights[ 3]   = Real(9.1490021622449999464462094123840e-02L);
        _weights[ 4]   = Real(1.1156664554733399471602390168177e-01L);
        _weights[ 5]   = Real(1.2875396253933622767551578485688e-01L);
        _weights[ 6]   = Real(1.4260670217360661177574610944190e-01L);
        _weights[ 7]   = Real(1.5276604206585966677885540089766e-01L);
        _weights[ 8]   = Real(1.5896884339395434764995643946505e-01L);
        _weights[ 9]   = Real(1.6105444984878369597916362532092e-01L);
        _weights[10]   = _weights[8];
        _weights[11]   = _weights[7];
        _weights[12]   = _weights[6];
        _weights[13]   = _weights[5];
        _weights[14]   = _weights[4];
        _weights[15]   = _weights[3];
        _weights[16]   = _weights[2];
        _weights[17]   = _weights[1];
        _weights[18]   = _weights[0];

        return;
      }

    case THIRTYEIGHTH:
    case THIRTYNINTH:
      {
        _points.resize (20);
        _weights.resize(20);

        _points[ 0](0) = Real(-9.9312859918509492478612238847132e-01L);
        _points[ 1](0) = Real(-9.6397192727791379126766613119728e-01L);
        _points[ 2](0) = Real(-9.1223442825132590586775244120330e-01L);
        _points[ 3](0) = Real(-8.3911697182221882339452906170152e-01L);
        _points[ 4](0) = Real(-7.4633190646015079261430507035564e-01L);
        _points[ 5](0) = Real(-6.3605368072651502545283669622629e-01L);
        _points[ 6](0) = Real(-5.1086700195082709800436405095525e-01L);
        _points[ 7](0) = Real(-3.7370608871541956067254817702493e-01L);
        _points[ 8](0) = Real(-2.2778585114164507808049619536857e-01L);
        _points[ 9](0) = Real(-7.6526521133497333754640409398838e-02L);
        _points[10]    = -_points[9];
        _points[11]    = -_points[8];
        _points[12]    = -_points[7];
        _points[13]    = -_points[6];
        _points[14]    = -_points[5];
        _points[15]    = -_points[4];
        _points[16]    = -_points[3];
        _points[17]    = -_points[2];
        _points[18]    = -_points[1];
        _points[19]    = -_points[0];

        _weights[ 0]   = Real(1.7614007139152118311861962351853e-02L);
        _weights[ 1]   = Real(4.0601429800386941331039952274932e-02L);
        _weights[ 2]   = Real(6.2672048334109063569506535187042e-02L);
        _weights[ 3]   = Real(8.3276741576704748724758143222046e-02L);
        _weights[ 4]   = Real(1.0193011981724043503675013548035e-01L);
        _weights[ 5]   = Real(1.1819453196151841731237737771138e-01L);
        _weights[ 6]   = Real(1.3168863844917662689849449974816e-01L);
        _weights[ 7]   = Real(1.4209610931838205132929832506716e-01L);
        _weights[ 8]   = Real(1.4917298647260374678782873700197e-01L);
        _weights[ 9]   = Real(1.5275338713072585069808433195510e-01L);
        _weights[10]   = _weights[9];
        _weights[11]   = _weights[8];
        _weights[12]   = _weights[7];
        _weights[13]   = _weights[6];
        _weights[14]   = _weights[5];
        _weights[15]   = _weights[4];
        _weights[16]   = _weights[3];
        _weights[17]   = _weights[2];
        _weights[18]   = _weights[1];
        _weights[19]   = _weights[0];

        return;
      }

    case FORTIETH:
    case FORTYFIRST:
      {
        _points.resize (21);
        _weights.resize(21);

        _points[ 0](0) = Real(-9.9375217062038950026024203593794e-01L);
        _points[ 1](0) = Real(-9.6722683856630629431662221490770e-01L);
        _points[ 2](0) = Real(-9.2009933415040082879018713371497e-01L);
        _points[ 3](0) = Real(-8.5336336458331728364725063858757e-01L);
        _points[ 4](0) = Real(-7.6843996347567790861587785130623e-01L);
        _points[ 5](0) = Real(-6.6713880419741231930596666999034e-01L);
        _points[ 6](0) = Real(-5.5161883588721980705901879672431e-01L);
        _points[ 7](0) = Real(-4.2434212020743878357366888854379e-01L);
        _points[ 8](0) = Real(-2.8802131680240109660079251606460e-01L);
        _points[ 9](0) = Real(-1.4556185416089509093703098233869e-01L);
        _points[10](0) = 0.;
        _points[11]    = -_points[9];
        _points[12]    = -_points[8];
        _points[13]    = -_points[7];
        _points[14]    = -_points[6];
        _points[15]    = -_points[5];
        _points[16]    = -_points[4];
        _points[17]    = -_points[3];
        _points[18]    = -_points[2];
        _points[19]    = -_points[1];
        _points[20]    = -_points[0];

        _weights[ 0]   = Real(1.6017228257774333324224616858471e-02L);
        _weights[ 1]   = Real(3.6953789770852493799950668299330e-02L);
        _weights[ 2]   = Real(5.7134425426857208283635826472448e-02L);
        _weights[ 3]   = Real(7.6100113628379302017051653300183e-02L);
        _weights[ 4]   = Real(9.3444423456033861553289741113932e-02L);
        _weights[ 5]   = Real(1.0879729916714837766347457807011e-01L);
        _weights[ 6]   = Real(1.2183141605372853419536717712572e-01L);
        _weights[ 7]   = Real(1.3226893863333746178105257449678e-01L);
        _weights[ 8]   = Real(1.3988739479107315472213342386758e-01L);
        _weights[ 9]   = Real(1.4452440398997005906382716655375e-01L);
        _weights[10]   = Real(1.4608113364969042719198514768337e-01L);
        _weights[11]   = _weights[9];
        _weights[12]   = _weights[8];
        _weights[13]   = _weights[7];
        _weights[14]   = _weights[6];
        _weights[15]   = _weights[5];
        _weights[16]   = _weights[4];
        _weights[17]   = _weights[3];
        _weights[18]   = _weights[2];
        _weights[19]   = _weights[1];
        _weights[20]   = _weights[0];

        return;
      }

    case FORTYSECOND:
    case FORTYTHIRD:
      {
        _points.resize (22);
        _weights.resize(22);

        _points[ 0](0) = Real(-9.9429458548239929207303142116130e-01L);
        _points[ 1](0) = Real(-9.7006049783542872712395098676527e-01L);
        _points[ 2](0) = Real(-9.2695677218717400052069293925905e-01L);
        _points[ 3](0) = Real(-8.6581257772030013653642563701938e-01L);
        _points[ 4](0) = Real(-7.8781680597920816200427795540835e-01L);
        _points[ 5](0) = Real(-6.9448726318668278005068983576226e-01L);
        _points[ 6](0) = Real(-5.8764040350691159295887692763865e-01L);
        _points[ 7](0) = Real(-4.6935583798675702640633071096641e-01L);
        _points[ 8](0) = Real(-3.4193582089208422515814742042738e-01L);
        _points[ 9](0) = Real(-2.0786042668822128547884653391955e-01L);
        _points[10](0) = Real(-6.9739273319722221213841796118628e-02L);
        _points[11]    = -_points[10];
        _points[12]    = -_points[9];
        _points[13]    = -_points[8];
        _points[14]    = -_points[7];
        _points[15]    = -_points[6];
        _points[16]    = -_points[5];
        _points[17]    = -_points[4];
        _points[18]    = -_points[3];
        _points[19]    = -_points[2];
        _points[20]    = -_points[1];
        _points[21]    = -_points[0];

        _weights[ 0]   = Real(1.4627995298272200684991098047185e-02L);
        _weights[ 1]   = Real(3.3774901584814154793302246865913e-02L);
        _weights[ 2]   = Real(5.2293335152683285940312051273211e-02L);
        _weights[ 3]   = Real(6.9796468424520488094961418930218e-02L);
        _weights[ 4]   = Real(8.5941606217067727414443681372703e-02L);
        _weights[ 5]   = Real(1.0041414444288096493207883783054e-01L);
        _weights[ 6]   = Real(1.1293229608053921839340060742175e-01L);
        _weights[ 7]   = Real(1.2325237681051242428556098615481e-01L);
        _weights[ 8]   = Real(1.3117350478706237073296499253031e-01L);
        _weights[ 9]   = Real(1.3654149834601517135257383123152e-01L);
        _weights[10]   = Real(1.3925187285563199337541024834181e-01L);
        _weights[11]   = _weights[10];
        _weights[12]   = _weights[9];
        _weights[13]   = _weights[8];
        _weights[14]   = _weights[7];
        _weights[15]   = _weights[6];
        _weights[16]   = _weights[5];
        _weights[17]   = _weights[4];
        _weights[18]   = _weights[3];
        _weights[19]   = _weights[2];
        _weights[20]   = _weights[1];
        _weights[21]   = _weights[0];

        return;
      }


    default:
      libmesh_error_msg("Quadrature rule " << _order << " not supported!");
    }
}

} // namespace libMesh
