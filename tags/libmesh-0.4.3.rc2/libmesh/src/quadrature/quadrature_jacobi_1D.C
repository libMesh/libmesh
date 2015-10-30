// $Id: quadrature_jacobi_1D.C,v 1.5 2004-01-03 15:37:44 benkirk Exp $

// The libMesh Finite Element Library.
// Copyright (C) 2002-2004  Benjamin S. Kirk, John W. Peterson
  
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



// Local includes
#include "quadrature_jacobi.h"
 


void QJacobi::init_1D(const ElemType)
{
  //----------------------------------------------------------------------
  // 1D quadrature rules
  // Note that we need to check _alpha and _beta
  // to determine which quadrature rule to implement.
  // These rules have been pre-scaled based on their
  // intended use.  The weights of the (alpha=1,beta=0)
  // rule sum to 0.5.  The weights of the (alpha=2,beta=0)
  // rule sum to 0.333333333333.

  if ((_alpha == 1) && (_beta == 0))
    {
      switch(_order)
	{
	case CONSTANT:
	case FIRST:
	  {
	    _points.resize (1);
	    _weights.resize(1);
	
	    _points[0](0)  = 0.33333333333333; 
	
	    _weights[0]    = 0.50000000000000;
	
	    return;
	  }
	case SECOND:
	case THIRD:
	  {
	    _points.resize (2);
	    _weights.resize(2);
	      
	    _points[0](0) = 0.15505102572168;
	    _points[1](0) = 0.64494897427832;

	    _weights[0]   = 0.31804138174398;
	    _weights[1]   = 0.18195861825602;

	    return;
	  }
	case FOURTH:
	case FIFTH:
	  {
	    _points.resize (3);
	    _weights.resize(3);
	      
	    _points[0](0) = 0.08858795951270;
	    _points[1](0) = 0.40946686444073;
	    _points[2](0) = 0.78765946176085;

	    _weights[0]   = 0.20093191373896;
	    _weights[1]   = 0.22924110635959;
	    _weights[2]   = 0.06982697990145;

	    return;
	  }
	case SIXTH:
	case SEVENTH:
	  {
	    _points.resize (4);
	    _weights.resize(4);
	      
	    _points[0](0) = 0.05710419611452;
	    _points[1](0) = 0.27684301363812;
	    _points[2](0) = 0.58359043236892;
	    _points[3](0) = 0.86024013565622;

	    _weights[0]   = 0.13550691343149;
	    _weights[1]   = 0.20346456801027;
	    _weights[2]   = 0.12984754760823;
	    _weights[3]   = 0.03118097095001;

	    return;
	  }
	case EIGHTH:
	case NINTH:
	  {
	    _points.resize (5);
	    _weights.resize(5);
	      
	    _points[0](0) = 0.03980985705147;
	    _points[1](0) = 0.19801341787361;
	    _points[2](0) = 0.43797481024739;
	    _points[3](0) = 0.69546427335364;
	    _points[4](0) = 0.90146491420117;
	      
	    _weights[0]   = 0.09678159022665;
	    _weights[1]   = 0.16717463809437;
	    _weights[2]   = 0.14638698708467;
	    _weights[3]   = 0.07390887007262;
	    _weights[4]   = 0.01574791452169;

	    return;
	  }
	case TENTH:
	case ELEVENTH:
	  {
	
	    _points.resize (6);
	    _weights.resize(6);
	      
	    _points[0](0) = 0.02931642715979;
	    _points[1](0) = 0.14807859966848;
	    _points[2](0) = 0.33698469028115;
	    _points[3](0) = 0.55867151877155;
	    _points[4](0) = 0.76923386203005;
	    _points[5](0) = 0.92694567131974;

	    _weights[0]   = 0.07231033072551;
	    _weights[1]   = 0.13554249723152;
	    _weights[2]   = 0.14079255378820;
	    _weights[3]   = 0.09866115089066;
	    _weights[4]   = 0.04395516555051;
	    _weights[5]   = 0.00873830181361;

	    return;
	  }
	case TWELFTH:
	case THIRTEENTH:
	  {
	    _points.resize (7);
	    _weights.resize(7);
	      
	    _points[0](0) = 0.02247938643871;
	    _points[1](0) = 0.11467905316090;
	    _points[2](0) = 0.26578982278459;
	    _points[3](0) = 0.45284637366944;
	    _points[4](0) = 0.64737528288683;
	    _points[5](0) = 0.81975930826311;
	    _points[6](0) = 0.94373743946308;

	    _weights[0]   = 0.05596736342349;
	    _weights[1]   = 0.11050925819087;
	    _weights[2]   = 0.12739089729959;
	    _weights[3]   = 0.10712506569587;
	    _weights[4]   = 0.06638469646549;
	    _weights[5]   = 0.02740835672187;
	    _weights[6]   = 0.00521436220281;

	    return;
	  }
	case FOURTEENTH:
	case FIFTEENTH:
	  {
	    _points.resize (8);
	    _weights.resize(8);
	      
	    _points[0](0) = 0.01777991514737;
	    _points[1](0) = 0.09132360789979;
	    _points[2](0) = 0.21430847939563;
	    _points[3](0) = 0.37193216458327;
	    _points[4](0) = 0.54518668480343;
	    _points[5](0) = 0.71317524285557;
	    _points[6](0) = 0.85563374295785;
	    _points[7](0) = 0.95536604471003;

	    _weights[0]   = 0.04455080436154;
	    _weights[1]   = 0.09111902363638;
	    _weights[2]   = 0.11250579947089;
	    _weights[3]   = 0.10604735943593;
	    _weights[4]   = 0.07919959949232;
	    _weights[5]   = 0.04543931950470;
	    _weights[6]   = 0.01784290265599;
	    _weights[7]   = 0.00329519144225;

	    return;
	  }
	case SIXTEENTH:
	case SEVENTEENTH:
	  {
	    _points.resize (9);
	    _weights.resize(9);

	    _points[0](0) = 0.01441240964887;
	    _points[1](0) = 0.07438738970920;
	    _points[2](0) = 0.17611665616299;
	    _points[3](0) = 0.30966757992764;
	    _points[4](0) = 0.46197040108101;
	    _points[5](0) = 0.61811723469529;
	    _points[6](0) = 0.76282301518504;
	    _points[7](0) = 0.88192102121000;
	    _points[8](0) = 0.96374218711679;

	    _weights[0]   = 0.03627800352333;
	    _weights[1]   = 0.07607425510930;
	    _weights[2]   = 0.09853374217235;
	    _weights[3]   = 0.10030880919337;
	    _weights[4]   = 0.08435832184492;
	    _weights[5]   = 0.05840119529517;
	    _weights[6]   = 0.03180482149105;
	    _weights[7]   = 0.01206000428479;
	    _weights[8]   = 0.00218084708577;

	    return;
	  }
	case EIGHTTEENTH:
	case NINTEENTH:
	  {
	    _points.resize (10);
	    _weights.resize(10);

	    _points[0](0) = 0.01191761343242;
	    _points[1](0) = 0.06173207187714;
	    _points[2](0) = 0.14711144964308;
	    _points[3](0) = 0.26115967600846;
	    _points[4](0) = 0.39463984688579;
	    _points[5](0) = 0.53673876571566;
	    _points[6](0) = 0.67594446167666;
	    _points[7](0) = 0.80097892103690;
	    _points[8](0) = 0.90171098779015;
	    _points[9](0) = 0.96997096783851;

	    _weights[0]   = 0.03009950802395;
	    _weights[1]   = 0.06428715450909;
	    _weights[2]   = 0.08621130028917;
	    _weights[3]   = 0.09269689367772;
	    _weights[4]   = 0.08455710969083;
	    _weights[5]   = 0.06605307556335;
	    _weights[6]   = 0.04340190640715;
	    _weights[7]   = 0.02277459145326;
	    _weights[8]   = 0.00841931978298;
	    _weights[9]   = 0.00149914060241;

	    return;
	  }      
	case TWENTIETH:
	case TWENTYFIRST:
	  {
	    _points.resize (11);
	    _weights.resize(11);

	    _points[0](0)  = 0.01001828046168;
	    _points[1](0)  = 0.05203545112718;
	    _points[2](0)  = 0.12461922514445;
	    _points[3](0)  = 0.22284060704384;
	    _points[4](0)  = 0.34000815791467;
	    _points[5](0)  = 0.46813761308958;
	    _points[6](0)  = 0.59849727976714;
	    _points[7](0)  = 0.72220328489097;
	    _points[8](0)  = 0.83082489962282;
	    _points[9](0)  = 0.91695838655260;
	    _points[10](0) = 0.97472637960248;

	    _weights[0]   = 0.02536734068817;
	    _weights[1]   = 0.05493809113287;
	    _weights[2]   = 0.07562004805718;
	    _weights[3]   = 0.08465942288402;
	    _weights[4]   = 0.08187910298806;
	    _weights[5]   = 0.06953187515818;
	    _weights[6]   = 0.05159136067230;
	    _weights[7]   = 0.03264154671383;
	    _weights[8]   = 0.01666362345168;
	    _weights[9]   = 0.00604392096048;
	    _weights[10]  = 0.00106366729324;

	    return;
	  }
	case TWENTYSECOND:
	case TWENTYTHIRD:
	  {
	    _points.resize (12);
	    _weights.resize(12);

	    _points[0](0)  = 0.00853905498844;
	    _points[1](0)  = 0.04444646315539;
	    _points[2](0)  = 0.10685449088348;
	    _points[3](0)  = 0.19215105452985;
	    _points[4](0)  = 0.29538088426258;
	    _points[5](0)  = 0.41054508120146;
	    _points[6](0)  = 0.53095084931282;
	    _points[7](0)  = 0.64960065027725;
	    _points[8](0)  = 0.75959888952523;
	    _points[9](0)  = 0.85455254376493;
	    _points[10](0) = 0.92894210126442;
	    _points[11](0) = 0.97843793683415;

	    _weights[0]   = 0.02166486088692;
	    _weights[1]   = 0.04742785198044;
	    _weights[2]   = 0.06660675062670;
	    _weights[3]   = 0.07689660268004;
	    _weights[4]   = 0.07769631681553;
	    _weights[5]   = 0.07010933999763;
	    _weights[6]   = 0.05661384371367;
	    _weights[7]   = 0.04045165370691;
	    _weights[8]   = 0.02487678040927;
	    _weights[9]   = 0.01243600916642;
	    _weights[10]  = 0.00444480779567;
	    _weights[11]  = 0.00077518222094;

	    return;
	  }      
	default:
	  {
	    std::cerr << "Quadrature rule " << _order
		      << " not supported!" << std::endl;
	      
	    error();
	  }
	}

      error();
    }


  
  else if ((_alpha == 2) && (_beta == 0))
    {
      
      switch(_order)
	{
	case CONSTANT:
	case FIRST:
	  {
	    _points.resize (1);
	    _weights.resize(1);
	
	    _points[0](0)  = 0.25000000000000;
	
	    _weights[0]    = 0.33333333333333;
	
	    return;
	  }
	case SECOND:
	case THIRD:
	  {
	    _points.resize (2);
	    _weights.resize(2);
	      
	    _points[0](0) = 0.12251482265544;
	    _points[1](0) = 0.54415184401123;

	    _weights[0]   = 0.23254745125351;
	    _weights[1]   = 0.10078588207983;

	    return;
	  }
	case FOURTH:
	case FIFTH:
	  {
	    _points.resize (3);
	    _weights.resize(3);
	      
	    _points[0](0) = 0.07299402407315;
	    _points[1](0) = 0.34700376603835;
	    _points[2](0) = 0.70500220988850;

	    _weights[0]   = 0.15713636106489;
	    _weights[1]   = 0.14624626925987;
	    _weights[2]   = 0.02995070300858;

	    return;
	  }
	case SIXTH:
	case SEVENTH:
	  {
	    _points.resize (4);
	    _weights.resize(4);
	      
	    _points[0](0) = 0.04850054944700;
	    _points[1](0) = 0.23860073755186;
	    _points[2](0) = 0.51704729510437;
	    _points[3](0) = 0.79585141789677;

	    _weights[0]   = 0.11088841561128;
	    _weights[1]   = 0.14345878979921;
	    _weights[2]   = 0.06863388717292;
	    _weights[3]   = 0.01035224074992;

	    return;
	  }
	case EIGHTH:
	case NINTH:
	  {
	    _points.resize (5);
	    _weights.resize(5);
	      
	    _points[0](0) = 0.03457893991821;
	    _points[1](0) = 0.17348032077170;
	    _points[2](0) = 0.38988638706552;
	    _points[3](0) = 0.63433347263089;
	    _points[4](0) = 0.85105421294702;
	      
	    _weights[0]   = 0.08176478428577;
	    _weights[1]   = 0.12619896189991;
	    _weights[2]   = 0.08920016122159;
	    _weights[3]   = 0.03205560072296;
	    _weights[4]   = 0.00411382520310;

	    return;
	  }
	case TENTH:
	case ELEVENTH:
	  {
	
	    _points.resize (6);
	    _weights.resize(6);
	      
	    _points[0](0) = 0.02590455509367;
	    _points[1](0) = 0.13156394165798;
	    _points[2](0) = 0.30243691802289;
	    _points[3](0) = 0.50903641316475;
	    _points[4](0) = 0.71568112731171;
	    _points[5](0) = 0.88680561617756;

	    _weights[0]   = 0.06253870272658;
	    _weights[1]   = 0.10737649973678;
	    _weights[2]   = 0.09457718674854;
	    _weights[3]   = 0.05128957112962;
	    _weights[4]   = 0.01572029718495;
	    _weights[5]   = 0.00183107580687;

	    return;
	  }
	case TWELFTH:
	case THIRTEENTH:
	  {
	    _points.resize (7);
	    _weights.resize(7);
	      
	    _points[0](0) = 0.02013277377340;
	    _points[1](0) = 0.10308902914805;
	    _points[2](0) = 0.24055412604806;
	    _points[3](0) = 0.41400214459706;
	    _points[4](0) = 0.60002151327899;
	    _points[5](0) = 0.77351724659144;
	    _points[6](0) = 0.91118316656300;

	    _weights[0]   = 0.04927650177644;
	    _weights[1]   = 0.09069882461269;
	    _weights[2]   = 0.09173380327980;
	    _weights[3]   = 0.06314637870889;
	    _weights[4]   = 0.02942221128953;
	    _weights[5]   = 0.00816292563230;
	    _weights[6]   = 0.00089268803369;

	    return;
	  }
	case FOURTEENTH:
	case FIFTEENTH:
	  {
	    _points.resize (8);
	    _weights.resize(8);
	      
	    _points[0](0) = 0.01609775955192;
	    _points[1](0) = 0.08290061748565;
	    _points[2](0) = 0.19547516848874;
	    _points[3](0) = 0.34165199147720;
	    _points[4](0) = 0.50559707818449;
	    _points[5](0) = 0.66955227182436;
	    _points[6](0) = 0.81577170358328;
	    _points[7](0) = 0.92850896495991;

	    _weights[0]   = 0.03977895780670;
	    _weights[1]   = 0.07681809326722;
	    _weights[2]   = 0.08528476917194;
	    _weights[3]   = 0.06844718342165;
	    _weights[4]   = 0.04081442638854;
	    _weights[5]   = 0.01724686378023;
	    _weights[6]   = 0.00447452171301;
	    _weights[7]   = 0.00046851778403;

	    return;
	  }
	case SIXTEENTH:
	case SEVENTEENTH:
	  {
	    _points.resize (9);
	    _weights.resize(9);

	    _points[0](0) = 0.01316588559711;
	    _points[1](0) = 0.06808452959377;
	    _points[2](0) = 0.16175951676407;
	    _points[3](0) = 0.28589108833922;
	    _points[4](0) = 0.42945364538781;
	    _points[5](0) = 0.57969405635116;
	    _points[6](0) = 0.72326857174034;
	    _points[7](0) = 0.84743684201324;
	    _points[8](0) = 0.94124586421327;

	    _weights[0]   = 0.03276014511105;
	    _weights[1]   = 0.06548953703338;
	    _weights[2]   = 0.07767356916056;
	    _weights[3]   = 0.06928439568980;
	    _weights[4]   = 0.04854062786451;
	    _weights[5]   = 0.02634328090255;
	    _weights[6]   = 0.01040611657935;
	    _weights[7]   = 0.00257439864561;
	    _weights[8]   = 0.00026126234652;

	    return;
	  }
	case EIGHTTEENTH:
	case NINTEENTH:
	  {
	    _points.resize (10);
	    _weights.resize(10);

	    _points[0](0) = 0.01096845245617;
	    _points[1](0) = 0.05689815053366;
	    _points[2](0) = 0.13595023405023;
	    _points[3](0) = 0.24228119613252;
	    _points[4](0) = 0.36800785044934;
	    _points[5](0) = 0.50380712641487;
	    _points[6](0) = 0.63960948865471;
	    _points[7](0) = 0.76534767954811;
	    _points[8](0) = 0.87171007457441;
	    _points[9](0) = 0.95087429264052;

	    _weights[0]   = 0.02743408871016;
	    _weights[1]   = 0.05627293640278;
	    _weights[2]   = 0.07006950770867;
	    _weights[3]   = 0.06745221938144;
	    _weights[4]   = 0.05288378876696;
	    _weights[5]   = 0.03385456501681;
	    _weights[6]   = 0.01719757504655;
	    _weights[7]   = 0.00646988906856;
	    _weights[8]   = 0.00154552319474;
	    _weights[9]   = 0.00015324003670;

	    return;
	  }      
	case TWENTIETH:
	case TWENTYFIRST:
	  {
	    _points.resize (11);
	    _weights.resize(11);

	    _points[0](0)  = 0.00927897383134;
	    _points[1](0)  = 0.04824969209430;
	    _points[2](0)  = 0.11578862662939;
	    _points[3](0)  = 0.20766834159706;
	    _points[4](0)  = 0.31811795190623;
	    _points[5](0)  = 0.44019839985886;
	    _points[6](0)  = 0.56623983915457;
	    _points[7](0)  = 0.68832423986296;
	    _points[8](0)  = 0.79878435859091;
	    _points[9](0)  = 0.89069109935439;
	    _points[10](0) = 0.95832514378665;

	    _weights[0]   = 0.02330085005155;
	    _weights[1]   = 0.04874586103051;
	    _weights[2]   = 0.06297954300041;
	    _weights[3]   = 0.06418355934975;
	    _weights[4]   = 0.05463222232105;
	    _weights[5]   = 0.03929021779844;
	    _weights[6]   = 0.02358966353881;
	    _weights[7]   = 0.01141459091810;
	    _weights[8]   = 0.00414000426322;
	    _weights[9]   = 0.00096302057554;
	    _weights[10]  = 0.00009380048601;

	    return;
	  }
	case TWENTYSECOND:
	case TWENTYTHIRD:
	  {
	    _points.resize (12);
	    _weights.resize(12);

	    _points[0](0)  = 0.00795204570266;
	    _points[1](0)  = 0.04142781045426;
	    _points[2](0)  = 0.09975762554264;
	    _points[3](0)  = 0.17981078905241;
	    _points[4](0)  = 0.27727345779932;
	    _points[5](0)  = 0.38689200999769;
	    _points[6](0)  = 0.50275736044903;
	    _points[7](0)  = 0.61862386345846;
	    _points[8](0)  = 0.72824645295307;
	    _points[9](0)  = 0.82571851421479;
	    _points[10](0) = 0.90579507354454;
	    _points[11](0) = 0.96420653529267;

	    _weights[0]   = 0.02003111258445;
	    _weights[1]   = 0.04255661536742;
	    _weights[2]   = 0.05658418474943;
	    _weights[3]   = 0.06025060748762;
	    _weights[4]   = 0.05457394116423;
	    _weights[5]   = 0.04276476357019;
	    _weights[6]   = 0.02890805183873;
	    _weights[7]   = 0.01654734364461;
	    _weights[8]   = 0.00771636702389;
	    _weights[9]   = 0.00272084214410;
	    _weights[10]  = 0.00061995339854;
	    _weights[11]  = 0.00005955035981;

	    return;
	  }      
	default:
	  {
	    std::cerr << "Quadrature rule " << _order
		      << " not supported!" << std::endl;
	      
	    error();
	  }
	}

      error();
    }

  else
    {
      std::cerr << "Unsupported combination of (alpha,beta) = ("
		<< _alpha
		<< ","
		<< _beta
		<< ") requested in Jacobi-Gauss quadrature rule."
		<< std::endl;
    }

  
  return;
}
 

