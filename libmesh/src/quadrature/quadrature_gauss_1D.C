// $Id: quadrature_gauss_1D.C,v 1.14 2006-03-22 19:26:33 roystgnr Exp $

// The libMesh Finite Element Library.
// Copyright (C) 2002-2005  Benjamin S. Kirk, John W. Peterson
  
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
#include "quadrature_gauss.h"
 


void QGauss::init_1D(const ElemType)
{
  //----------------------------------------------------------------------
  // 1D quadrature rules
  switch(_order)
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
	      
	_points[0](0) = -0.577350269189626; // -sqrt(3)/3
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
	      
	_points[0](0) = -0.774596669241483377035853079956;
	_points[1](0) = 0.;
	_points[2]    = -_points[0];

	_weights[0]   = 0.555555555555555555555555555556;
	_weights[1]   = 0.888888888888888888888888888889;
	_weights[2]   = _weights[0];

	return;
      }
    case SIXTH:
    case SEVENTH:
      {
	_points.resize (4);
	_weights.resize(4);
	      
	_points[0](0) = -0.861136311594052575223946488893;
	_points[1](0) = -0.339981043584856264802665759103;
	_points[2]    = -_points[1];
	_points[3]    = -_points[0];

	_weights[0]   = 0.347854845137453857373063949222;
	_weights[1]   = 0.652145154862546142626936050778;
	_weights[2]   = _weights[1];
	_weights[3]   = _weights[0];

	return;
      }
    case EIGHTH:
    case NINTH:
      {
	_points.resize (5);
	_weights.resize(5);
	      
	_points[0](0) = -0.906179845938663992797626878299;
	_points[1](0) = -0.538469310105683091036314420700;
	_points[2](0) = 0.;
	_points[3]    = -_points[1];
	_points[4]    = -_points[0];
	      
	_weights[0]   = 0.236926885056189087514264040720;
	_weights[1]   = 0.478628670499366468041291514836;
	_weights[2]   = 0.568888888888888888888888888889;
	_weights[3]   = _weights[1];
	_weights[4]   = _weights[0];

	return;
      }
    case TENTH:
    case ELEVENTH:
      {
	_points.resize (6);
	_weights.resize(6);
	      
	_points[0](0) = -0.932469514203152027812301554494;
	_points[1](0) = -0.661209386466264513661399595020;
	_points[2](0) = -0.238619186083196908630501721681;
	_points[3]    = -_points[2];
	_points[4]    = -_points[1];
	_points[5]    = -_points[0];

	_weights[0]   = 0.171324492379170345040296142173;
	_weights[1]   = 0.360761573048138607569833513838;
	_weights[2]   = 0.467913934572691047389870343990;
	_weights[3]   = _weights[2];
	_weights[4]   = _weights[1];
	_weights[5]   = _weights[0];

	return;
      }
    case TWELFTH:
    case THIRTEENTH:
      {
	_points.resize (7);
	_weights.resize(7);
	      
	_points[0](0) = -0.949107912342758524526189684048;
	_points[1](0) = -0.741531185599394439863864773281;
	_points[2](0) = -0.405845151377397166906606412077;
	_points[3](0) = 0.;
	_points[4]    = -_points[2];
	_points[5]    = -_points[1];
	_points[6]    = -_points[0];

	_weights[0]   = 0.12948496616887;
	_weights[1]   = 0.27970539148928;
	_weights[2]   = 0.38183005050512;
	_weights[3]   = 0.41795918367347;
	_weights[4]   = _weights[2];
	_weights[5]   = _weights[1];
	_weights[6]   = _weights[0];

	return;
      }
    case FOURTEENTH:
    case FIFTEENTH:
      {
	_points.resize (8);
	_weights.resize(8);
	      
	_points[0](0) = -0.960289856497536231683560868569;
	_points[1](0) = -0.796666477413626739591553936476;
	_points[2](0) = -0.525532409916328985817739049189;
	_points[3](0) = -0.183434642495649804939476142360;
	_points[4]    = -_points[3];
	_points[5]    = -_points[2];
	_points[6]    = -_points[1];
	_points[7]    = -_points[0];

	_weights[0]   = 0.101228536290376259152531354310;
	_weights[1]   = 0.222381034453374470544355994426;
	_weights[2]   = 0.313706645877887287337962201987;
	_weights[3]   = 0.362683783378361982965150449277;
	_weights[4]   = _weights[3];
	_weights[5]   = _weights[2];
	_weights[6]   = _weights[1];
	_weights[7]   = _weights[0];

	return;
      }
    case SIXTEENTH:
    case SEVENTEENTH:
      {
	_points.resize (9);
	_weights.resize(9);

	_points[0](0) = -0.968160239507626089835576202904;
	_points[1](0) = -0.836031107326635794299429788070;
	_points[2](0) = -0.613371432700590397308702039341;
	_points[3](0) = -0.324253423403808929038538014643;
	_points[4](0) =  0.0000000000000000000000000000000;
	_points[5]    = -_points[3];
	_points[6]    = -_points[2];
	_points[7]    = -_points[1];
	_points[8]    = -_points[0];

	_weights[0]   = 0.0812743883615744119718921581105;
	_weights[1]   = 0.180648160694857404058472031243;
	_weights[2]   = 0.260610696402935462318742869419;
	_weights[3]   = 0.312347077040002840068630406584;
	_weights[4]   = 0.330239355001259763164525069287;
	_weights[5]   = _weights[3];
	_weights[6]   = _weights[2];
	_weights[7]   = _weights[1];
	_weights[8]   = _weights[0];

	return;
      }
    case EIGHTTEENTH:
    case NINTEENTH:
      {
	_points.resize (10);
	_weights.resize(10);

	_points[0](0) = -0.973906528517171720077964012084;
	_points[1](0) = -0.865063366688984510732096688423;
	_points[2](0) = -0.679409568299024406234327365115;
	_points[3](0) = -0.433395394129247190799265943166;
	_points[4](0) = -0.148874338981631210864826001130;
	_points[5]    = -_points[4];
	_points[6]    = -_points[3];
	_points[7]    = -_points[2];
	_points[8]    = -_points[1];
	_points[9]    = -_points[0];

	_weights[0]   = 0.0666713443086881375935688098933; 
	_weights[1]   = 0.149451349150580593145776339658;
	_weights[2]   = 0.219086362515982043995534934228;
	_weights[3]   = 0.269266719309996355091226921569;
	_weights[4]   = 0.295524224714752870173892994651;
	_weights[5]   = _weights[4];
	_weights[6]   = _weights[3];
	_weights[7]   = _weights[2];
	_weights[8]   = _weights[1];
	_weights[9]   = _weights[0];

	return;
      }      

      /*
       * Note that from Order TWENTIETH upwards
       * (except TWENTYSECOND and TWENTYTHIRD),
       * the 14th digit may not always be correct,
       * only _close_ to the 14th digit of the
       * true coordinate and weight...
       */

    case TWENTIETH:
    case TWENTYFIRST:
      {
	_points.resize (11);
	_weights.resize(11);

	_points[0](0)  = -0.97822865814606;
	_points[1](0)  = -0.88706259976810;
	_points[2](0)  = -0.73015200557405;
	_points[3](0)  = -0.51909612920681;
	_points[4](0)  = -0.26954315595235;
	_points[5](0)  = 0.;
	_points[6]     = -_points[4];
	_points[7]     = -_points[3];
	_points[8]     = -_points[2];
	_points[9]     = -_points[1];
	_points[10]    = -_points[0];

	_weights[0]     = 0.05566856711617;
	_weights[1]     = 0.12558036946490;
	_weights[2]     = 0.18629021092773;
	_weights[3]     = 0.23319376459199;
	_weights[4]     = 0.26280454451025;
	_weights[5]     = 0.27292508677790;
	_weights[6]     = _weights[4];
	_weights[7]     = _weights[3];
	_weights[8]     = _weights[2];
	_weights[9]     = _weights[1];
	_weights[10]    = _weights[0];

	return;
      }

    case TWENTYSECOND:
    case TWENTYTHIRD:
      {
	_points.resize (12);
	_weights.resize(12);

	_points[0](0) = -0.981560634246719250690549090149;
	_points[1](0) = -0.904117256370474856678465866119;
	_points[2](0) = -0.769902674194304687036893833213;
	_points[3](0) = -0.587317954286617447296702418941;
	_points[4](0) = -0.367831498998180193752691536644;
	_points[5](0) = -0.125233408511468915472441369464;
	_points[6]    = -_points[5];
	_points[7]    = -_points[4];
	_points[8]    = -_points[3];
	_points[9]    = -_points[2];
	_points[10]   = -_points[1];
	_points[11]   = -_points[0];

	_weights[0]   = 0.0471753363865118271946159614850; 
	_weights[1]   = 0.106939325995318430960254718194;
	_weights[2]   = 0.160078328543346226334652529543;
	_weights[3]   = 0.203167426723065921749064455810;
	_weights[4]   = 0.233492536538354808760849898925;
	_weights[5]   = 0.249147045813402785000562436043;
	_weights[6]   = _weights[5];
	_weights[7]   = _weights[4];
	_weights[8]   = _weights[3];
	_weights[9]   = _weights[2];
	_weights[10]  = _weights[1];
	_weights[11]  = _weights[0];

	return;
      }      

    case TWENTYFOURTH:
    case TWENTYFIFTH:
      {
	_points.resize (13);
	_weights.resize(13);

	_points[0](0)  = -0.98418305471859;
	_points[1](0)  = -0.91759839922298;
	_points[2](0)  = -0.80157809073331;
	_points[3](0)  = -0.64234933944034;
	_points[4](0)  = -0.44849275103645;
	_points[5](0)  = -0.23045831595513;
	_points[6](0)  = 0.;
	_points[7]     = -_points[5];
	_points[8]     = -_points[4];
	_points[9]     = -_points[3];
	_points[10]    = -_points[2];
	_points[11]    = -_points[1];
	_points[12]    = -_points[0];

	_weights[0]     = 0.04048400476532;
	_weights[1]     = 0.09212149983773;
	_weights[2]     = 0.13887351021979;
	_weights[3]     = 0.17814598076195;
	_weights[4]     = 0.20781604753689;
	_weights[5]     = 0.22628318026290;
	_weights[6]     = 0.23255155323087;
	_weights[7]     = _weights[5];
	_weights[8]     = _weights[4];
	_weights[9]     = _weights[3];
	_weights[10]    = _weights[2];
	_weights[11]    = _weights[1];
	_weights[12]    = _weights[0];

	return;
      }

    case TWENTYSIXTH:
    case TWENTYSEVENTH:
      {
	_points.resize (14);
	_weights.resize(14);

	_points[0](0)  = -0.98628380869681;
	_points[1](0)  = -0.92843488366357;
	_points[2](0)  = -0.82720131506977;
	_points[3](0)  = -0.68729290481169;
	_points[4](0)  = -0.51524863635815;
	_points[5](0)  = -0.31911236892789;
	_points[6](0)  = -0.10805494870734;
	_points[7]     = -_points[6];
	_points[8]     = -_points[5];
	_points[9]     = -_points[4];
	_points[10]    = -_points[3];
	_points[11]    = -_points[2];
	_points[12]    = -_points[1];
	_points[13]    = -_points[0];

	_weights[0]     = 0.03511946033175;
	_weights[1]     = 0.08015808715976;
	_weights[2]     = 0.12151857068790;
	_weights[3]     = 0.15720316715819;
	_weights[4]     = 0.18553839747794;
	_weights[5]     = 0.20519846372130;
	_weights[6]     = 0.21526385346316;
	_weights[7]     = _weights[6];
	_weights[8]     = _weights[5];
	_weights[9]     = _weights[4];
	_weights[10]    = _weights[3];
	_weights[11]    = _weights[2];
	_weights[12]    = _weights[1];
	_weights[13]    = _weights[0];

	return;
      }

    case TWENTYEIGHTH:
    case TWENTYNINTH:
      {
	_points.resize (15);
	_weights.resize(15);

	_points[0](0)  = -0.98799251802049;
	_points[1](0)  = -0.93727339240071;
	_points[2](0)  = -0.84820658341043;
	_points[3](0)  = -0.72441773136017;
	_points[4](0)  = -0.57097217260854;
	_points[5](0)  = -0.39415134707756;
	_points[6](0)  = -0.20119409399743;
	_points[7](0)  = -0.00000000000000;
	_points[8]     = -_points[6];
	_points[9]     = -_points[5];
	_points[10]    = -_points[4];
	_points[11]    = -_points[3];
	_points[12]    = -_points[2];
	_points[13]    = -_points[1];
	_points[14]    = -_points[0];

	_weights[0]     = 0.03075324199612;
	_weights[1]     = 0.07036604748811;
	_weights[2]     = 0.10715922046717;
	_weights[3]     = 0.13957067792615;
	_weights[4]     = 0.16626920581699;
	_weights[5]     = 0.18616100001556;
	_weights[6]     = 0.19843148532711;
	_weights[7]     = 0.20257824192556;
	_weights[8]     = _weights[6];
	_weights[9]     = _weights[5];
	_weights[10]    = _weights[4];
	_weights[11]    = _weights[3];
	_weights[12]    = _weights[2];
	_weights[13]    = _weights[1];
	_weights[14]    = _weights[0];

	return;
      }

    case THIRTIETH:
    case THIRTYFIRST:
      {
	_points.resize (16);
	_weights.resize(16);

	_points[0](0)  = -0.98940093499165;
	_points[1](0)  = -0.94457502307323;
	_points[2](0)  = -0.86563120238783;
	_points[3](0)  = -0.75540440835500;
	_points[4](0)  = -0.61787624440264;
	_points[5](0)  = -0.45801677765723;
	_points[6](0)  = -0.28160355077926;
	_points[7](0)  = -0.09501250983764;
	_points[8]     = -_points[7];
	_points[9]     = -_points[6];
	_points[10]    = -_points[5];
	_points[11]    = -_points[4];
	_points[12]    = -_points[3];
	_points[13]    = -_points[2];
	_points[14]    = -_points[1];
	_points[15]    = -_points[0];

	_weights[0]     = 0.02715245941175;
	_weights[1]     = 0.06225352393865;
	_weights[2]     = 0.09515851168249;
	_weights[3]     = 0.12462897125553;
	_weights[4]     = 0.14959598881658;
	_weights[5]     = 0.16915651939500;
	_weights[6]     = 0.18260341504492;
	_weights[7]     = 0.18945061045507;
	_weights[8]     = _weights[7];
	_weights[9]     = _weights[6];
	_weights[10]    = _weights[5];
	_weights[11]    = _weights[4];
	_weights[12]    = _weights[3];
	_weights[13]    = _weights[2];
	_weights[14]    = _weights[1];
	_weights[15]    = _weights[0];

	return;
      }

    case THIRTYSECOND:
    case THIRTYTHIRD:
      {
	_points.resize (17);
	_weights.resize(17);

	_points[0](0)  = -0.99057547531442;
	_points[1](0)  = -0.95067552176877;
	_points[2](0)  = -0.88023915372699;
	_points[3](0)  = -0.78151400389680;
	_points[4](0)  = -0.65767115921669;
	_points[5](0)  = -0.51269053708648;
	_points[6](0)  = -0.35123176345388;
	_points[7](0)  = -0.17848418149585;
	_points[8](0)  = 0.;
	_points[9]     = -_points[7];
	_points[10]    = -_points[6];
	_points[11]    = -_points[5];
	_points[12]    = -_points[4];
	_points[13]    = -_points[3];
	_points[14]    = -_points[2];
	_points[15]    = -_points[1];
	_points[16]    = -_points[0];

	_weights[0]     = 0.02414830286855;
	_weights[1]     = 0.05545952937399;
	_weights[2]     = 0.08503614831718;
	_weights[3]     = 0.11188384719340;
	_weights[4]     = 0.13513636846853;
	_weights[5]     = 0.15404576107681;
	_weights[6]     = 0.16800410215645;
	_weights[7]     = 0.17656270536699;
	_weights[8]     = 0.17944647035621;
	_weights[9]     = _weights[7];
	_weights[10]    = _weights[6];
	_weights[11]    = _weights[5];
	_weights[12]    = _weights[4];
	_weights[13]    = _weights[3];
	_weights[14]    = _weights[2];
	_weights[15]    = _weights[1];
	_weights[16]    = _weights[0];

	return;
      }

    case THIRTYFOURTH:
    case THIRTYFIFTH:
      {
	_points.resize (18);
	_weights.resize(18);

	_points[0](0)  = -0.99156516842093;
	_points[1](0)  = -0.95582394957140;
	_points[2](0)  = -0.89260246649756;
	_points[3](0)  = -0.80370495897252;
	_points[4](0)  = -0.69168704306035;
	_points[5](0)  = -0.55977083107395;
	_points[6](0)  = -0.41175116146284;
	_points[7](0)  = -0.25188622569151;
	_points[8](0)  = -0.08477501304173;
	_points[9]     = -_points[8];
	_points[10]    = -_points[7];
	_points[11]    = -_points[6];
	_points[12]    = -_points[5];
	_points[13]    = -_points[4];
	_points[14]    = -_points[3];
	_points[15]    = -_points[2];
	_points[16]    = -_points[1];
	_points[17]    = -_points[0];

	_weights[0]     = 0.02161601352648;
	_weights[1]     = 0.04971454889497;
	_weights[2]     = 0.07642573025489;
	_weights[3]     = 0.10094204410629;
	_weights[4]     = 0.12255520671148;
	_weights[5]     = 0.14064291467065;
	_weights[6]     = 0.15468467512627;
	_weights[7]     = 0.16427648374583;
	_weights[8]     = 0.16914238296314;
	_weights[9]     = _weights[8];
	_weights[10]    = _weights[7];
	_weights[11]    = _weights[6];
	_weights[12]    = _weights[5];
	_weights[13]    = _weights[4];
	_weights[14]    = _weights[3];
	_weights[15]    = _weights[2];
	_weights[16]    = _weights[1];
	_weights[17]    = _weights[0];

	return;
      }

    case THIRTYSIXTH:
    case THIRTYSEVENTH:
      {
	_points.resize (19);
	_weights.resize(19);

	_points[0](0)  = -0.99240684384358;
	_points[1](0)  = -0.96020815213483;
	_points[2](0)  = -0.90315590361482;
	_points[3](0)  = -0.82271465653714;
	_points[4](0)  = -0.72096617733523;
	_points[5](0)  = -0.60054530466168;
	_points[6](0)  = -0.46457074137596;
	_points[7](0)  = -0.31656409996363;
	_points[8](0)  = -0.16035864564023;
	_points[9](0)  = 0.;
	_points[10]    = -_points[8];
	_points[11]    = -_points[7];
	_points[12]    = -_points[6];
	_points[13]    = -_points[5];
	_points[14]    = -_points[4];
	_points[15]    = -_points[3];
	_points[16]    = -_points[2];
	_points[17]    = -_points[1];
	_points[18]    = -_points[0];

	_weights[0]     = 0.01946178822973;
	_weights[1]     = 0.04481422676570;
	_weights[2]     = 0.06904454273764;
	_weights[3]     = 0.09149002162245;
	_weights[4]     = 0.11156664554733;
	_weights[5]     = 0.12875396253934;
	_weights[6]     = 0.14260670217361;
	_weights[7]     = 0.15276604206586;
	_weights[8]     = 0.15896884339395;
	_weights[9]     = 0.16105444984878;
	_weights[10]    = _weights[8];
	_weights[11]    = _weights[7];
	_weights[12]    = _weights[6];
	_weights[13]    = _weights[5];
	_weights[14]    = _weights[4];
	_weights[15]    = _weights[3];
	_weights[16]    = _weights[2];
	_weights[17]    = _weights[1];
	_weights[18]    = _weights[0];

	return;
      }

    case THIRTYEIGHTH:
    case THIRTYNINTH:
      {
	_points.resize (20);
	_weights.resize(20);

	_points[0](0)  = -0.99312859918510;
	_points[1](0)  = -0.96397192727791;
	_points[2](0)  = -0.91223442825133;
	_points[3](0)  = -0.83911697182222;
	_points[4](0)  = -0.74633190646015;
	_points[5](0)  = -0.63605368072652;
	_points[6](0)  = -0.51086700195083;
	_points[7](0)  = -0.37370608871542;
	_points[8](0)  = -0.22778585114164;
	_points[9](0)  = -0.07652652113350;
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

	_weights[0]     = 0.01761400713915;
	_weights[1]     = 0.04060142980039;
	_weights[2]     = 0.06267204833411;
	_weights[3]     = 0.08327674157670;
	_weights[4]     = 0.10193011981724;
	_weights[5]     = 0.11819453196152;
	_weights[6]     = 0.13168863844918;
	_weights[7]     = 0.14209610931838;
	_weights[8]     = 0.14917298647260;
	_weights[9]     = 0.15275338713073;
	_weights[10]    = _weights[9];
	_weights[11]    = _weights[8];
	_weights[12]    = _weights[7];
	_weights[13]    = _weights[6];
	_weights[14]    = _weights[5];
	_weights[15]    = _weights[4];
	_weights[16]    = _weights[3];
	_weights[17]    = _weights[2];
	_weights[18]    = _weights[1];
	_weights[19]    = _weights[0];

	return;
      }

    case FORTIETH:
    case FORTYFIRST:
      {
	_points.resize (21);
	_weights.resize(21);

	_points[0](0)  = -0.99375217062039;
	_points[1](0)  = -0.96722683856631;
	_points[2](0)  = -0.92009933415040;
	_points[3](0)  = -0.85336336458332;
	_points[4](0)  = -0.76843996347568;
	_points[5](0)  = -0.66713880419741;
	_points[6](0)  = -0.55161883588722;
	_points[7](0)  = -0.42434212020744;
	_points[8](0)  = -0.28802131680240;
	_points[9](0)  = -0.14556185416090;
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

	_weights[0]     = 0.01601722825777;
	_weights[1]     = 0.03695378977085;
	_weights[2]     = 0.05713442542686;
	_weights[3]     = 0.07610011362838;
	_weights[4]     = 0.09344442345603;
	_weights[5]     = 0.10879729916715;
	_weights[6]     = 0.12183141605373;
	_weights[7]     = 0.13226893863334;
	_weights[8]     = 0.13988739479107;
	_weights[9]     = 0.14452440398997;
	_weights[10]    = 0.14608113364969;
	_weights[11]    = _weights[9];
	_weights[12]    = _weights[8];
	_weights[13]    = _weights[7];
	_weights[14]    = _weights[6];
	_weights[15]    = _weights[5];
	_weights[16]    = _weights[4];
	_weights[17]    = _weights[3];
	_weights[18]    = _weights[2];
	_weights[19]    = _weights[1];
	_weights[20]    = _weights[0];

	return;
      }

    case FORTYSECOND:
    case FORTYTHIRD:
      {
	_points.resize (22);
	_weights.resize(22);

	_points[0](0)  = -0.99429458548240;
	_points[1](0)  = -0.97006049783543;
	_points[2](0)  = -0.92695677218717;
	_points[3](0)  = -0.86581257772030;
	_points[4](0)  = -0.78781680597921;
	_points[5](0)  = -0.69448726318668;
	_points[6](0)  = -0.58764040350691;
	_points[7](0)  = -0.46935583798676;
	_points[8](0)  = -0.34193582089208;
	_points[9](0)  = -0.20786042668822;
	_points[10](0) = -0.06973927331972;
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

	_weights[0]     = 0.01462799529827;
	_weights[1]     = 0.03377490158481;
	_weights[2]     = 0.05229333515268;
	_weights[3]     = 0.06979646842452;
	_weights[4]     = 0.08594160621707;
	_weights[5]     = 0.10041414444288;
	_weights[6]     = 0.11293229608054;
	_weights[7]     = 0.12325237681051;
	_weights[8]     = 0.13117350478706;
	_weights[9]     = 0.13654149834602;
	_weights[10]    = 0.13925187285563;
	_weights[11]    = _weights[10];
	_weights[12]    = _weights[9];
	_weights[13]    = _weights[8];
	_weights[14]    = _weights[7];
	_weights[15]    = _weights[6];
	_weights[16]    = _weights[5];
	_weights[17]    = _weights[4];
	_weights[18]    = _weights[3];
	_weights[19]    = _weights[2];
	_weights[20]    = _weights[1];
	_weights[21]    = _weights[0];

	return;
      }


    default:
      {
	std::cerr << "Quadrature rule " << _order
		  << " not supported!" << std::endl;
	      
	error();
      }
    }



  return;
}
 

