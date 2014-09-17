// The libMesh Finite Element Library.
// Copyright (C) 2002-2014 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

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
#include "libmesh/quadrature_gauss_lobatto.h"

namespace libMesh
{

void QGaussLobatto::init_1D(const ElemType,
                            unsigned int p)
{
  //----------------------------------------------------------------------
  // 1D quadrature rules
  switch(_order + 2*p)
    {
      // Since Gauss-Lobatto rules must include the endpoints of the
      // domain, there is no 1-point rule.  The two-point
      // Gauss-Lobatto rule is equivalent to the trapezoidal rule.
    case CONSTANT:
    case FIRST:
      {
        _points.resize (2);
        _weights.resize(2);

        _points[0](0) = -1.0L;
        _points[1]    = -_points[0];

        _weights[0]   = 1.;
        _weights[1]   = _weights[0];

        return;
      }

      // The three-point Gauss-Lobatto rule is equivalent to Simpsons' rule.
      // It can integrate cubic polynomials exactly.
    case SECOND:
    case THIRD:
      {
        _points.resize (3);
        _weights.resize(3);

        _points[0](0) = -1.0L;
        _points[1]    = 0.0L;
        _points[2]    = -_points[0];

        _weights[0]   = 1.0L / 3.0L;
        _weights[1]   = 4.0L / 3.0L;
        _weights[2]   = _weights[0];
        return;
      }

      // The four-point Gauss-Lobatto rule can integrate 2*4-3 =
      // 5th-order polynomials exactly.
    case FOURTH:
    case FIFTH:
      {
        _points.resize (4);
        _weights.resize(4);

        _points[ 0](0) = -1.0L;
        _points[ 1](0) = -std::sqrt(1.0L/5.0L);
        _points[ 2]    = -_points[1];
        _points[ 3]    = -_points[0];

        _weights[ 0]   = 1.0L/6.0L;
        _weights[ 1]   = 5.0L/6.0L;
        _weights[ 2]   = _weights[1];
        _weights[ 3]   = _weights[0];

        return;
      }

      // The five-point Gauss-Lobatto rule can integrate 2*5-3 =
      // 7th-order polynomials exactly.
    case SIXTH:
    case SEVENTH:
      {
        _points.resize (5);
        _weights.resize(5);

        _points[ 0](0) = -1.0L;
        _points[ 1](0) = -std::sqrt(3.0L/7.0L);
        _points[ 2](0) = 0.;
        _points[ 3]    = -_points[1];
        _points[ 4]    = -_points[0];

        _weights[ 0]   = 1.0L/10.0L;
        _weights[ 1]   = 49.0L/90.0L;
        _weights[ 2]   = 32.0L/45.0L;
        _weights[ 3]   = _weights[1];
        _weights[ 4]   = _weights[0];

        return;
      }

      // 2*6-3 = 9
    case EIGHTH:
    case NINTH:
      {
        _points.resize (6);
        _weights.resize(6);

        _points[ 0](0) = -1.0L;
        _points[ 1](0) = -7.65055323929464692851002973959e-01L;
        _points[ 2](0) = -2.85231516480645096314150994041e-01L;
        _points[ 3]    = -_points[2];
        _points[ 4]    = -_points[1];
        _points[ 5]    = -_points[0];

        _weights[ 0]   = 1.0L/15.0L;
        _weights[ 1]   = 3.7847495629784698031661280821e-01L;
        _weights[ 2]   = 5.5485837703548635301672052512e-01L;
        _weights[ 3]   = _weights[2];
        _weights[ 4]   = _weights[1];
        _weights[ 5]   = _weights[0];

        return;
      }

      // 2*7-3 = 11
    case TENTH:
    case ELEVENTH:
      {
        _points.resize (7);
        _weights.resize(7);

        _points[ 0](0) = -1.0L;
        _points[ 1](0) = -8.30223896278566929872032213967e-01L;
        _points[ 2](0) = -4.68848793470714213803771881909e-01L;
        _points[ 3](0) = 0.;
        _points[ 4]    = -_points[2];
        _points[ 5]    = -_points[1];
        _points[ 6]    = -_points[0];

        _weights[ 0]   = 1.0L/21.0L;
        _weights[ 1]   = 2.7682604736156594801070040629e-01L;
        _weights[ 2]   = 4.3174538120986262341787102228e-01L;
        _weights[ 3]   = 4.8761904761904761904761904762e-01L;
        _weights[ 4]   = _weights[2];
        _weights[ 5]   = _weights[1];
        _weights[ 6]   = _weights[0];

        return;
      }

      // 2*8-3 = 13
    case TWELFTH:
    case THIRTEENTH:
      {
        _points.resize (8);
        _weights.resize(8);

        _points[ 0](0) = -1.0L;
        _points[ 1](0) = -8.71740148509606615337445761221e-01L;
        _points[ 2](0) = -5.91700181433142302144510731398e-01L;
        _points[ 3](0) = -2.09299217902478868768657260345e-01L;
        _points[ 4]    = -_points[3];
        _points[ 5]    = -_points[2];
        _points[ 6]    = -_points[1];
        _points[ 7]    = -_points[0];

        _weights[ 0]   = 1.0L/28.0L;
        _weights[ 1]   = 2.107042271435060393829920658e-01L;
        _weights[ 2]   = 3.411226924835043647642406771e-01L;
        _weights[ 3]   = 4.124587946587038815670529714e-01L;
        _weights[ 4]   = _weights[3];
        _weights[ 5]   = _weights[2];
        _weights[ 6]   = _weights[1];
        _weights[ 7]   = _weights[0];

        return;
      }

      // 2*9-3 = 15
    case FOURTEENTH:
    case FIFTEENTH:
      {
        _points.resize (9);
        _weights.resize(9);

        _points[ 0](0) = -1.0L;
        _points[ 1](0) = -8.99757995411460157312345244418e-01L;
        _points[ 2](0) = -6.77186279510737753445885427091e-01L;
        _points[ 3](0) = -3.63117463826178158710752068709e-01L;
        _points[ 4](0) = 0.;
        _points[ 5]    = -_points[3];
        _points[ 6]    = -_points[2];
        _points[ 7]    = -_points[1];
        _points[ 8]    = -_points[0];

        _weights[ 0]   = 1.0L/36.0L;
        _weights[ 1]   = 1.6549536156080552504633972003e-01L;
        _weights[ 2]   = 2.7453871250016173528070561858e-01L;
        _weights[ 3]   = 3.4642851097304634511513153214e-01L;
        _weights[ 4]   = 3.7151927437641723356009070295e-01L;
        _weights[ 5]   = _weights[3];
        _weights[ 6]   = _weights[2];
        _weights[ 7]   = _weights[1];
        _weights[ 8]   = _weights[0];

        return;
      }

      // 2*10-3 = 17
    case SIXTEENTH:
    case SEVENTEENTH:
      {
        _points.resize (10);
        _weights.resize(10);

        _points[ 0](0) = -1.0L;
        _points[ 1](0) = -9.1953390816645881382893266082e-01L;
        _points[ 2](0) = -7.3877386510550507500310617486e-01L;
        _points[ 3](0) = -4.7792494981044449566117509273e-01L;
        _points[ 4](0) = -1.6527895766638702462621976596e-01L;
        _points[ 5]    = -_points[4];
        _points[ 6]    = -_points[3];
        _points[ 7]    = -_points[2];
        _points[ 8]    = -_points[1];
        _points[ 9]    = -_points[0];

        _weights[ 0]   = 1.0L/45.0L;
        _weights[ 1]   = 1.3330599085107011112622717076e-01L;
        _weights[ 2]   = 2.2488934206312645211945782173e-01L;
        _weights[ 3]   = 2.9204268367968375787558225737e-01L;
        _weights[ 4]   = 3.2753976118389745665651052792e-01L;
        _weights[ 5]   = _weights[4];
        _weights[ 6]   = _weights[3];
        _weights[ 7]   = _weights[2];
        _weights[ 8]   = _weights[1];
        _weights[ 9]   = _weights[0];

        return;
      }

      // 2*11-3 = 19
    case EIGHTTEENTH:
    case NINETEENTH:
      {
        _points.resize (11);
        _weights.resize(11);

        _points[ 0](0) = -1.0L;
        _points[ 1](0) = -9.34001430408059134332274136099e-01L;
        _points[ 2](0) = -7.84483473663144418622417816108e-01L;
        _points[ 3](0) = -5.65235326996205006470963969478e-01L;
        _points[ 4](0) = -2.95758135586939391431911515559e-01L;
        _points[ 5](0) = 0.;
        _points[ 6]    = -_points[4];
        _points[ 7]    = -_points[3];
        _points[ 8]    = -_points[2];
        _points[ 9]    = -_points[1];
        _points[10]    = -_points[0];

        _weights[ 0]   = 1.0L/55.0L;
        _weights[ 1]   = 1.096122732669948644614034496e-01L;
        _weights[ 2]   = 1.871698817803052041081415219e-01L;
        _weights[ 3]   = 2.480481042640283140400848664e-01L;
        _weights[ 4]   = 2.868791247790080886792224033e-01L;
        _weights[ 5]   = 3.002175954556906937859318812e-01L;
        _weights[ 6]   = _weights[4];
        _weights[ 7]   = _weights[3];
        _weights[ 8]   = _weights[2];
        _weights[ 9]   = _weights[1];
        _weights[10]   = _weights[0];

        return;
      }

      // 2*12-3 = 21
    case TWENTIETH:
    case TWENTYFIRST:
      {
        _points.resize (12);
        _weights.resize(12);

        _points[ 0](0) = -1.0L;
        _points[ 1](0) = -9.4489927222288222340758013830e-01L;
        _points[ 2](0) = -8.1927932164400667834864158172e-01L;
        _points[ 3](0) = -6.3287615303186067766240485444e-01L;
        _points[ 4](0) = -3.9953094096534893226434979157e-01L;
        _points[ 5](0) = -1.3655293285492755486406185574e-01L;
        _points[ 6]    = -_points[5];
        _points[ 7]    = -_points[4];
        _points[ 8]    = -_points[3];
        _points[ 9]    = -_points[2];
        _points[10]    = -_points[1];
        _points[11]    = -_points[0];

        _weights[ 0]   = 1.0L/66.0L;
        _weights[ 1]   = 9.168451741319613066834259413e-02L;
        _weights[ 2]   = 1.579747055643701151646710627e-01L;
        _weights[ 3]   = 2.125084177610211453583020774e-01L;
        _weights[ 4]   = 2.512756031992012802932444122e-01L;
        _weights[ 5]   = 2.714052409106961770002883385e-01L;
        _weights[ 6]   = _weights[5];
        _weights[ 7]   = _weights[4];
        _weights[ 8]   = _weights[3];
        _weights[ 9]   = _weights[2];
        _weights[10]   = _weights[1];
        _weights[11]   = _weights[0];

        return;
      }

      // 2*13-3 = 23
    case TWENTYSECOND:
    case TWENTYTHIRD:
      {
        _points.resize (13);
        _weights.resize(13);

        _points[ 0](0) = -1.0L;
        _points[ 1](0) = -9.53309846642163911896905464755e-01L;
        _points[ 2](0) = -8.46347564651872316865925607099e-01L;
        _points[ 3](0) = -6.86188469081757426072759039566e-01L;
        _points[ 4](0) = -4.82909821091336201746937233637e-01L;
        _points[ 5](0) = -2.49286930106239992568673700374e-01L;
        _points[ 6](0) = 0.;
        _points[ 7]    = -_points[5];
        _points[ 8]    = -_points[4];
        _points[ 9]    = -_points[3];
        _points[10]    = -_points[2];
        _points[11]    = -_points[1];
        _points[12]    = -_points[0];

        _weights[ 0]   = 1.0/78.0L;
        _weights[ 1]   = 7.7801686746818927793588988333e-02L;
        _weights[ 2]   = 1.3498192668960834911991476259e-01L;
        _weights[ 3]   = 1.8364686520355009200749425875e-01L;
        _weights[ 4]   = 2.2076779356611008608553400838e-01L;
        _weights[ 5]   = 2.4401579030667635645857814836e-01L;
        _weights[ 6]   = 2.5193084933344673604413864154e-01L;
        _weights[ 7]   = _weights[5];
        _weights[ 8]   = _weights[4];
        _weights[ 9]   = _weights[3];
        _weights[10]   = _weights[2];
        _weights[11]   = _weights[1];
        _weights[12]   = _weights[0];

        return;
      }

      // 2*14-3 = 25
    case TWENTYFOURTH:
    case TWENTYFIFTH:
      {
        _points.resize (14);
        _weights.resize(14);

        _points[ 0](0) = -1.0L;
        _points[ 1](0) = -9.5993504526726090135510016202e-01L;
        _points[ 2](0) = -8.6780105383034725100022020291e-01L;
        _points[ 3](0) = -7.2886859909132614058467240052e-01L;
        _points[ 4](0) = -5.5063940292864705531662270586e-01L;
        _points[ 5](0) = -3.4272401334271284504390340364e-01L;
        _points[ 6](0) = -1.1633186888370386765877670974e-01L;
        _points[ 7]    = -_points[6];
        _points[ 8]    = -_points[5];
        _points[ 9]    = -_points[4];
        _points[10]    = -_points[3];
        _points[11]    = -_points[2];
        _points[12]    = -_points[1];
        _points[13]    = -_points[0];

        _weights[ 0]   = 1.0L/91.0L;
        _weights[ 1]   = 6.6837284497681284634070660746e-02L;
        _weights[ 2]   = 1.1658665589871165154099667066e-01L;
        _weights[ 3]   = 1.6002185176295214241282099799e-01L;
        _weights[ 4]   = 1.9482614937341611864033177838e-01L;
        _weights[ 5]   = 2.1912625300977075487116252395e-01L;
        _weights[ 6]   = 2.3161279446845705888962835729e-01L;
        _weights[ 7]   = _weights[6];
        _weights[ 8]   = _weights[5];
        _weights[ 9]   = _weights[4];
        _weights[10]   = _weights[3];
        _weights[11]   = _weights[2];
        _weights[12]   = _weights[1];
        _weights[13]   = _weights[0];

        return;
      }

      // 2*15-3 = 27
    case TWENTYSIXTH:
    case TWENTYSEVENTH:
      {
        _points.resize (15);
        _weights.resize(15);

        _points[ 0](0) = -1.0L;
        _points[ 1](0) = -9.6524592650383857279585139207e-01L;
        _points[ 2](0) = -8.8508204422297629882540163148e-01L;
        _points[ 3](0) = -7.6351968995181520070411847598e-01L;
        _points[ 4](0) = -6.0625320546984571112352993864e-01L;
        _points[ 5](0) = -4.2063805471367248092189693874e-01L;
        _points[ 6](0) = -2.1535395536379423822567944627e-01L;
        _points[ 7](0) = 0.;
        _points[ 8]    = -_points[6];
        _points[ 9]    = -_points[5];
        _points[10]    = -_points[4];
        _points[11]    = -_points[3];
        _points[12]    = -_points[2];
        _points[13]    = -_points[1];
        _points[14]    = -_points[0];

        _weights[ 0]   = 1.0L/105.L;
        _weights[ 1]   = 5.8029893028601249096880584025e-02L;
        _weights[ 2]   = 1.0166007032571806760366617079e-01L;
        _weights[ 3]   = 1.4051169980242810946044680564e-01L;
        _weights[ 4]   = 1.7278964725360094905207709941e-01L;
        _weights[ 5]   = 1.9698723596461335609250034651e-01L;
        _weights[ 6]   = 2.1197358592682092012743007698e-01L;
        _weights[ 7]   = 2.1704811634881564951495021425e-01L;
        _weights[ 8]   = _weights[6];
        _weights[ 9]   = _weights[5];
        _weights[10]   = _weights[4];
        _weights[11]   = _weights[3];
        _weights[12]   = _weights[2];
        _weights[13]   = _weights[1];
        _weights[14]   = _weights[0];

        return;
      }

      // 2*16-3 = 29
    case TWENTYEIGHTH:
    case TWENTYNINTH:
      {
        _points.resize (16);
        _weights.resize(16);

        _points[ 0](0) = -1.0L;
        _points[ 1](0) = -9.6956804627021793295224273837e-01L;
        _points[ 2](0) = -8.9920053309347209299462826152e-01L;
        _points[ 3](0) = -7.9200829186181506393108827096e-01L;
        _points[ 4](0) = -6.5238870288249308946788321964e-01L;
        _points[ 5](0) = -4.8605942188713761178189078585e-01L;
        _points[ 6](0) = -2.9983046890076320809835345472e-01L;
        _points[ 7](0) = -1.0132627352194944784303300505e-01L;
        _points[ 8]    = -_points[7];
        _points[ 9]    = -_points[6];
        _points[10]    = -_points[5];
        _points[11]    = -_points[4];
        _points[12]    = -_points[3];
        _points[13]    = -_points[2];
        _points[14]    = -_points[1];
        _points[15]    = -_points[0];

        _weights[ 0]   = 1.0L/120.0L;
        _weights[ 1]   = 5.0850361005919905403244919565e-02L;
        _weights[ 2]   = 8.9393697325930800991052080166e-02L;
        _weights[ 3]   = 1.2425538213251409834953633266e-01L;
        _weights[ 4]   = 1.5402698080716428081564494049e-01L;
        _weights[ 5]   = 1.7749191339170412530107566953e-01L;
        _weights[ 6]   = 1.9369002382520358431691359885e-01L;
        _weights[ 7]   = 2.0195830817822987148919912541e-01L;
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

      // 2*17-3 = 31
    case THIRTIETH:
    case THIRTYFIRST:
      {
        _points.resize (17);
        _weights.resize(17);

        _points[ 0](0) = -1.0L;
        _points[ 1](0) = -9.73132176631418314156979501874e-01L;
        _points[ 2](0) = -9.10879995915573595623802506398e-01L;
        _points[ 3](0) = -8.15696251221770307106750553238e-01L;
        _points[ 4](0) = -6.91028980627684705394919357372e-01L;
        _points[ 5](0) = -5.41385399330101539123733407504e-01L;
        _points[ 6](0) = -3.72174433565477041907234680735e-01L;
        _points[ 7](0) = -1.89511973518317388304263014753e-01L;
        _points[ 8](0) = 0.;
        _points[ 9]    = -_points[7];
        _points[10]    = -_points[6];
        _points[11]    = -_points[5];
        _points[12]    = -_points[4];
        _points[13]    = -_points[3];
        _points[14]    = -_points[2];
        _points[15]    = -_points[1];
        _points[16]    = -_points[0];

        _weights[ 0]   = 1.0L/136.0L;
        _weights[ 1]   = 4.4921940543254209647400954623e-02L;
        _weights[ 2]   = 7.9198270503687119190264429953e-02L;
        _weights[ 3]   = 1.1059290900702816137577270522e-01L;
        _weights[ 4]   = 1.3798774620192655905620157495e-01L;
        _weights[ 5]   = 1.6039466199762153951632836586e-01L;
        _weights[ 6]   = 1.7700425351565787043694574536e-01L;
        _weights[ 7]   = 1.8721633967761923589208848286e-01L;
        _weights[ 8]   = 1.9066187475346943329940724703e-01L;
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

      // 2*18-3 = 33
    case THIRTYSECOND:
    case THIRTYTHIRD:
      {
        _points.resize (18);
        _weights.resize(18);

        _points[ 0](0) = -1.0L;
        _points[ 1](0) = -9.76105557412198542864518924342e-01L;
        _points[ 2](0) = -9.20649185347533873837854625431e-01L;
        _points[ 3](0) = -8.35593535218090213713646362328e-01L;
        _points[ 4](0) = -7.23679329283242681306210365302e-01L;
        _points[ 5](0) = -5.88504834318661761173535893194e-01L;
        _points[ 6](0) = -4.34415036912123975342287136741e-01L;
        _points[ 7](0) = -2.66362652878280984167665332026e-01L;
        _points[ 8](0) = -8.97490934846521110226450100886e-02L;
        _points[ 9]    = -_points[8];
        _points[10]    = -_points[7];
        _points[11]    = -_points[6];
        _points[12]    = -_points[5];
        _points[13]    = -_points[4];
        _points[14]    = -_points[3];
        _points[15]    = -_points[2];
        _points[16]    = -_points[1];
        _points[17]    = -_points[0];

        _weights[ 0]   = 1.0L/153.0L;
        _weights[ 1]   = 3.997062881091406613759917641e-02L;
        _weights[ 2]   = 7.063716688563366499922296017e-02L;
        _weights[ 3]   = 9.901627171750280239442360532e-02L;
        _weights[ 4]   = 1.242105331329671002633963589e-01L;
        _weights[ 5]   = 1.454119615738022679830032105e-01L;
        _weights[ 6]   = 1.619395172376024892643267067e-01L;
        _weights[ 7]   = 1.732621094894562260106144038e-01L;
        _weights[ 8]   = 1.790158634397030822938188069e-01L;
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

      // 2*19-3 = 35
    case THIRTYFOURTH:
    case THIRTYFIFTH:
      {
        _points.resize (19);
        _weights.resize(19);

        _points[ 0](0) = -1.0L;
        _points[ 1](0) = -9.7861176622208009515263406311e-01L;
        _points[ 2](0) = -9.2890152815258624371794025880e-01L;
        _points[ 3](0) = -8.5246057779664609308595597004e-01L;
        _points[ 4](0) = -7.5149420255261301416363748963e-01L;
        _points[ 5](0) = -6.2890813726522049776683230623e-01L;
        _points[ 6](0) = -4.8822928568071350277790963763e-01L;
        _points[ 7](0) = -3.3350484782449861029850010385e-01L;
        _points[ 8](0) = -1.6918602340928157137515415345e-01L;
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

        _weights[ 0]   = 1.0L/171.0L;
        _weights[ 1]   = 3.579336518617647711542556904e-02L;
        _weights[ 2]   = 6.338189176262973685169569042e-02L;
        _weights[ 3]   = 8.913175709920708444800879056e-02L;
        _weights[ 4]   = 1.123153414773050440709100155e-01L;
        _weights[ 5]   = 1.322672804487507769260467339e-01L;
        _weights[ 6]   = 1.484139425959388850096806437e-01L;
        _weights[ 7]   = 1.602909240440612419799109682e-01L;
        _weights[ 8]   = 1.675565845271428672701372777e-01L;
        _weights[ 9]   = 1.700019192848272346446727156e-01L;
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

      // 2*20-3 = 37
    case THIRTYSIXTH:
    case THIRTYSEVENTH:
      {
        _points.resize (20);
        _weights.resize(20);

        _points[ 0](0) = -1.0L;
        _points[ 1](0) = -9.807437048939141719254464386e-01L;
        _points[ 2](0) = -9.359344988126654357161815849e-01L;
        _points[ 3](0) = -8.668779780899501413098472146e-01L;
        _points[ 4](0) = -7.753682609520558704143175276e-01L;
        _points[ 5](0) = -6.637764022903112898464033230e-01L;
        _points[ 6](0) = -5.349928640318862616481359618e-01L;
        _points[ 7](0) = -3.923531837139092993864747038e-01L;
        _points[ 8](0) = -2.395517059229864951824013569e-01L;
        _points[ 9](0) = -8.054593723882183797594451816e-02L;
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

        _weights[ 0]   = 1.0L/190.0L;
        _weights[ 1]   = 3.2237123188488941491605028117e-02L;
        _weights[ 2]   = 5.7181802127566826004753627173e-02L;
        _weights[ 3]   = 8.0631763996119603144776846114e-02L;
        _weights[ 4]   = 1.0199149969945081568378120573e-01L;
        _weights[ 5]   = 1.2070922762867472509942970500e-01L;
        _weights[ 6]   = 1.3630048235872418448978079299e-01L;
        _weights[ 7]   = 1.4836155407091682581471301373e-01L;
        _weights[ 8]   = 1.5658010264747548715816989679e-01L;
        _weights[ 9]   = 1.6074328638784574900772672645e-01L;
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

      // 2*21-3 = 39
    case THIRTYEIGHTH:
    case THIRTYNINTH:
      {
        _points.resize (21);
        _weights.resize(21);

        _points[ 0](0) = -1.0L;
        _points[ 1](0) = -9.8257229660454802823448127656e-01L;
        _points[ 2](0) = -9.4197629695974553429610265066e-01L;
        _points[ 3](0) = -8.7929475532359046445115359631e-01L;
        _points[ 4](0) = -7.9600192607771240474431258966e-01L;
        _points[ 5](0) = -6.9405102606222323262731639320e-01L;
        _points[ 6](0) = -5.7583196026183068692702187034e-01L;
        _points[ 7](0) = -4.4411578327900210119451634961e-01L;
        _points[ 8](0) = -3.0198985650876488727535186786e-01L;
        _points[ 9](0) = -1.5278551580218546600635832849e-01L;
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

        _weights[ 0]   = 1.0L/210.0L;
        _weights[ 1]   = 2.9184840098505458609458543613e-02L;
        _weights[ 2]   = 5.1843169000849625072722971853e-02L;
        _weights[ 3]   = 7.3273918185074144252547861042e-02L;
        _weights[ 4]   = 9.2985467957886065301137664149e-02L;
        _weights[ 5]   = 1.1051708321912333526700048678e-01L;
        _weights[ 6]   = 1.2545812119086894801515753571e-01L;
        _weights[ 7]   = 1.3745846286004134358089961742e-01L;
        _weights[ 8]   = 1.4623686244797745926727053063e-01L;
        _weights[ 9]   = 1.5158757511168138445325068151e-01L;
        _weights[10]   = 1.5338519033217494855158440507e-01L;
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

      // 2*22-3 = 41
    case FORTIETH:
    case FORTYFIRST:
      {
        _points.resize (22);
        _weights.resize(22);

        _points[ 0](0) = -1.0L;
        _points[ 1](0) = -9.8415243845764617655228962221e-01L;
        _points[ 2](0) = -9.4720428399922868052421376662e-01L;
        _points[ 3](0) = -8.9006229019090447052965782578e-01L;
        _points[ 4](0) = -8.1394892761192113604544184806e-01L;
        _points[ 5](0) = -7.2048723996120215811988189640e-01L;
        _points[ 6](0) = -6.1166943828425897122621160587e-01L;
        _points[ 7](0) = -4.8981487518990234980875123568e-01L;
        _points[ 8](0) = -3.5752071013891953806095728024e-01L;
        _points[ 9](0) = -2.1760658515928504178795509347e-01L;
        _points[10](0) = -7.3054540010898334761088790464e-02L;
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

        _weights[ 0]   = 1.0L/231.0L;
        _weights[ 1]   = 2.6545747682501757911627904521e-02L;
        _weights[ 2]   = 4.7214465293740752123775734865e-02L;
        _weights[ 3]   = 6.6865605864553076012404194157e-02L;
        _weights[ 4]   = 8.5090060391838447815711236096e-02L;
        _weights[ 5]   = 1.0150057480164767437243730375e-01L;
        _weights[ 6]   = 1.1574764465393906659003636772e-01L;
        _weights[ 7]   = 1.2752769665343027553084445931e-01L;
        _weights[ 8]   = 1.3658968861374142668617736221e-01L;
        _weights[ 9]   = 1.4274049227136140033623599357e-01L;
        _weights[10]   = 1.4584901944424179361642043948e-01L;
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

      // 2*23-3 = 43
    case FORTYSECOND:
    case FORTYTHIRD:
      {
        _points.resize (23);
        _weights.resize(23);

        _points[ 0](0) = -1.0L;
        _points[ 1](0) = -9.8552715587873257808146276674e-01L;
        _points[ 2](0) = -9.5175795571071020413563967985e-01L;
        _points[ 3](0) = -8.9945855804034501095016032035e-01L;
        _points[ 4](0) = -8.2965109665128588622320061929e-01L;
        _points[ 5](0) = -7.4369504117206068394516354307e-01L;
        _points[ 6](0) = -6.4326364446013620847614553360e-01L;
        _points[ 7](0) = -5.3031177113684416813011532015e-01L;
        _points[ 8](0) = -4.0703793791447482919595048822e-01L;
        _points[ 9](0) = -2.7584154894579306710687763268e-01L;
        _points[10](0) = -1.3927620404066839859186261298e-01L;
        _points[11]    = 0.;
        _points[12]    = -_points[10];
        _points[13]    = -_points[9];
        _points[14]    = -_points[8];
        _points[15]    = -_points[7];
        _points[16]    = -_points[6];
        _points[17]    = -_points[5];
        _points[18]    = -_points[4];
        _points[19]    = -_points[3];
        _points[20]    = -_points[2];
        _points[21]    = -_points[1];
        _points[22]    = -_points[0];

        _weights[ 0]   = 1.0L/253.0L;
        _weights[ 1]   = 2.424860077153173651739965894e-02L;
        _weights[ 2]   = 4.317587117024183474887646561e-02L;
        _weights[ 3]   = 6.125247712955420638138284744e-02L;
        _weights[ 4]   = 7.813544947556998974193425535e-02L;
        _weights[ 5]   = 9.349724616351234183350070691e-02L;
        _weights[ 6]   = 1.070391017243365115351836279e-01L;
        _weights[ 7]   = 1.184975106627491313021260047e-01L;
        _weights[ 8]   = 1.276494747017588766361485531e-01L;
        _weights[ 9]   = 1.343168726386038199015648977e-01L;
        _weights[10]   = 1.383699363858073945235027339e-01L;
        _weights[11]   = 1.397297800127473651401597065e-01L;
        _weights[12]   = _weights[10];
        _weights[13]   = _weights[9];
        _weights[14]   = _weights[8];
        _weights[15]   = _weights[7];
        _weights[16]   = _weights[6];
        _weights[17]   = _weights[5];
        _weights[18]   = _weights[4];
        _weights[19]   = _weights[3];
        _weights[20]   = _weights[2];
        _weights[21]   = _weights[1];
        _weights[22]   = _weights[0];

        return;
      }

    default:
      libmesh_error_msg("Quadrature rule " << _order << " not supported!");
    }
}

} // namespace libMesh
