// $Id: quadrature_gauss_1D.C,v 1.5 2003-02-06 06:02:42 jwpeterson Exp $

// The Next Great Finite Element Library.
// Copyright (C) 2002  Benjamin S. Kirk, John W. Peterson
  
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
    case CONST:
    case FIRST:
      {
	_points.resize (1);
	_weights.resize(1);
	
	_points[0](0)  = 0.;
	
	_weights[0]    = 2.;
	
	break;
      };
    case SECOND:
    case THIRD:
      {
	_points.resize (2);
	_weights.resize(2);
	      
	_points[0](0) = -0.577350269189626;
	_points[1]    = -_points[0];

	_weights[0]   = 1.;
	_weights[1]   = _weights[0];

	break;
      };
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

	break;
      };
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

	break;
      };
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

	break;
      };
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

	break;
      };
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

	break;
      };
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

	break;
      };
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

	break;
      };
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

	break;
      };      
    case TWENTIETH:
    case TWENTYFIRST:
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

	break;
      };      
    default:
      {
	std::cerr << "Quadrature rule " << _order
		  << " not supported!" << std::endl;
	      
	error();
      };
    };



  return;
};
 

