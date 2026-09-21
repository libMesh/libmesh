// The libMesh Finite Element Library.
// Copyright (C) 2002-2026 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

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
#include <array>
#include <vector>

// Local includes
#include "libmesh/int_range.h"
#include "libmesh/quadrature_gauss_lobatto.h"

namespace libMesh
{

namespace {

/**
 * Fills \p pts and \p wts with the Gauss-Lobatto points and weights of the rule of order
 * \p order on the reference interval [-1,1], the points in ascending order.
 *
 * A rule of n points is exact through degree 2n-3, so the two orders 2n-4 and 2n-3 share
 * a rule, which is why the cases below come in pairs.
 */
void fill_1D(const Order order,
             std::vector<Real> & pts,
             std::vector<Real> & wts)
{
  switch(order)
    {
      // Since Gauss-Lobatto rules must include the endpoints of the
      // domain, there is no 1-point rule.  The two-point
      // Gauss-Lobatto rule is equivalent to the trapezoidal rule.
    case CONSTANT:
    case FIRST:
      {
        pts.resize (2);
        wts.resize(2);

        pts[0] = -1;
        pts[1] = -pts[0];

        wts[0] = 1.;
        wts[1] = wts[0];

        return;
      }

      // The three-point Gauss-Lobatto rule is equivalent to Simpsons' rule.
      // It can integrate cubic polynomials exactly.
    case SECOND:
    case THIRD:
      {
        pts.resize (3);
        wts.resize(3);

        pts[0] = -1;
        pts[1] = 0;
        pts[2] = -pts[0];

        wts[0] = Real(1)/3;
        wts[1] = Real(4)/3;
        wts[2] = wts[0];
        return;
      }

      // The four-point Gauss-Lobatto rule can integrate 2*4-3 =
      // 5th-order polynomials exactly.
    case FOURTH:
    case FIFTH:
      {
        pts.resize (4);
        wts.resize(4);

        pts[ 0] = -1;
        pts[ 1] = -std::sqrt(Real(1)/5); // ~ -0.4472135955
        pts[ 2] = -pts[1];
        pts[ 3] = -pts[0];

        wts[ 0] = Real(1)/6;
        wts[ 1] = Real(5)/6;
        wts[ 2] = wts[1];
        wts[ 3] = wts[0];

        return;
      }

      // The five-point Gauss-Lobatto rule can integrate 2*5-3 =
      // 7th-order polynomials exactly.
    case SIXTH:
    case SEVENTH:
      {
        pts.resize (5);
        wts.resize(5);

        pts[ 0] = -1;
        pts[ 1] = -std::sqrt(Real(3)/7); // ~ -0.6546536707
        pts[ 2] = 0.;
        pts[ 3] = -pts[1];
        pts[ 4] = -pts[0];

        wts[ 0] = 1/Real(10);
        wts[ 1] = Real(49)/90;
        wts[ 2] = Real(32)/45;
        wts[ 3] = wts[1];
        wts[ 4] = wts[0];

        return;
      }

      // 2*6-3 = 9
    case EIGHTH:
    case NINTH:
      {
        pts.resize (6);
        wts.resize(6);

        pts[ 0] = -1.0000000000000000000000000000000e+00_R;
        pts[ 1] = -7.6505532392946469285100297395934e-01_R;
        pts[ 2] = -2.8523151648064509631415099404088e-01_R;
        pts[ 3] = -pts[2];
        pts[ 4] = -pts[1];
        pts[ 5] = -pts[0];

        wts[ 0] = 6.6666666666666666666666666666667e-02_R;
        wts[ 1] = 3.7847495629784698031661280821202e-01_R;
        wts[ 2] = 5.5485837703548635301672052512131e-01_R;
        wts[ 3] = wts[2];
        wts[ 4] = wts[1];
        wts[ 5] = wts[0];

        return;
      }

      // 2*7-3 = 11
    case TENTH:
    case ELEVENTH:
      {
        pts.resize (7);
        wts.resize(7);

        pts[ 0] = -1.0000000000000000000000000000000e+00_R;
        pts[ 1] = -8.3022389627856692987203221396747e-01_R;
        pts[ 2] = -4.6884879347071421380377188190877e-01_R;
        pts[ 3] = 0.;
        pts[ 4] = -pts[2];
        pts[ 5] = -pts[1];
        pts[ 6] = -pts[0];

        wts[ 0] = 4.7619047619047619047619047619048e-02_R;
        wts[ 1] = 2.7682604736156594801070040629007e-01_R;
        wts[ 2] = 4.3174538120986262341787102228136e-01_R;
        wts[ 3] = 4.8761904761904761904761904761905e-01_R;
        wts[ 4] = wts[2];
        wts[ 5] = wts[1];
        wts[ 6] = wts[0];

        return;
      }

      // 2*8-3 = 13
    case TWELFTH:
    case THIRTEENTH:
      {
        pts.resize (8);
        wts.resize(8);

        pts[ 0] = -1.0000000000000000000000000000000e+00_R;
        pts[ 1] = -8.7174014850960661533744576122066e-01_R;
        pts[ 2] = -5.9170018143314230214451073139795e-01_R;
        pts[ 3] = -2.0929921790247886876865726034535e-01_R;
        pts[ 4] = -pts[3];
        pts[ 5] = -pts[2];
        pts[ 6] = -pts[1];
        pts[ 7] = -pts[0];

        wts[ 0] = 3.5714285714285714285714285714286e-02_R;
        wts[ 1] = 2.1070422714350603938299206577576e-01_R;
        wts[ 2] = 3.4112269248350436476424067710775e-01_R;
        wts[ 3] = 4.1245879465870388156705297140221e-01_R;
        wts[ 4] = wts[3];
        wts[ 5] = wts[2];
        wts[ 6] = wts[1];
        wts[ 7] = wts[0];

        return;
      }

      // 2*9-3 = 15
    case FOURTEENTH:
    case FIFTEENTH:
      {
        pts.resize (9);
        wts.resize(9);

        pts[ 0] = -1.0000000000000000000000000000000e+00_R;
        pts[ 1] = -8.9975799541146015731234524441834e-01_R;
        pts[ 2] = -6.7718627951073775344588542709134e-01_R;
        pts[ 3] = -3.6311746382617815871075206870866e-01_R;
        pts[ 4] = 0.;
        pts[ 5] = -pts[3];
        pts[ 6] = -pts[2];
        pts[ 7] = -pts[1];
        pts[ 8] = -pts[0];

        wts[ 0] = 2.7777777777777777777777777777778e-02_R;
        wts[ 1] = 1.6549536156080552504633972002921e-01_R;
        wts[ 2] = 2.7453871250016173528070561857937e-01_R;
        wts[ 3] = 3.4642851097304634511513153213972e-01_R;
        wts[ 4] = 3.7151927437641723356009070294785e-01_R;
        wts[ 5] = wts[3];
        wts[ 6] = wts[2];
        wts[ 7] = wts[1];
        wts[ 8] = wts[0];

        return;
      }

      // 2*10-3 = 17
    case SIXTEENTH:
    case SEVENTEENTH:
      {
        pts.resize (10);
        wts.resize(10);

        pts[ 0] = -1.0000000000000000000000000000000e+00_R;
        pts[ 1] = -9.1953390816645881382893266082234e-01_R;
        pts[ 2] = -7.3877386510550507500310617485983e-01_R;
        pts[ 3] = -4.7792494981044449566117509273126e-01_R;
        pts[ 4] = -1.6527895766638702462621976595817e-01_R;
        pts[ 5] = -pts[4];
        pts[ 6] = -pts[3];
        pts[ 7] = -pts[2];
        pts[ 8] = -pts[1];
        pts[ 9] = -pts[0];

        wts[ 0] = 2.2222222222222222222222222222222e-02_R;
        wts[ 1] = 1.3330599085107011112622717075539e-01_R;
        wts[ 2] = 2.2488934206312645211945782173105e-01_R;
        wts[ 3] = 2.9204268367968375787558225737444e-01_R;
        wts[ 4] = 3.2753976118389745665651052791689e-01_R;
        wts[ 5] = wts[4];
        wts[ 6] = wts[3];
        wts[ 7] = wts[2];
        wts[ 8] = wts[1];
        wts[ 9] = wts[0];

        return;
      }

      // 2*11-3 = 19
    case EIGHTTEENTH:
    case NINETEENTH:
      {
        pts.resize (11);
        wts.resize(11);

        pts[ 0] = -1.0000000000000000000000000000000e+00_R;
        pts[ 1] = -9.3400143040805913433227413609938e-01_R;
        pts[ 2] = -7.8448347366314441862241781610846e-01_R;
        pts[ 3] = -5.6523532699620500647096396947775e-01_R;
        pts[ 4] = -2.9575813558693939143191151555906e-01_R;
        pts[ 5] = 0.;
        pts[ 6] = -pts[4];
        pts[ 7] = -pts[3];
        pts[ 8] = -pts[2];
        pts[ 9] = -pts[1];
        pts[10] = -pts[0];

        wts[ 0] = 1.8181818181818181818181818181818e-02_R;
        wts[ 1] = 1.0961227326699486446140344958035e-01_R;
        wts[ 2] = 1.8716988178030520410814152189943e-01_R;
        wts[ 3] = 2.4804810426402831404008486642187e-01_R;
        wts[ 4] = 2.8687912477900808867922240333154e-01_R;
        wts[ 5] = 3.0021759545569069378593188116998e-01_R;
        wts[ 6] = wts[4];
        wts[ 7] = wts[3];
        wts[ 8] = wts[2];
        wts[ 9] = wts[1];
        wts[10] = wts[0];

        return;
      }

      // 2*12-3 = 21
    case TWENTIETH:
    case TWENTYFIRST:
      {
        pts.resize (12);
        wts.resize(12);

        pts[ 0] = -1.0000000000000000000000000000000e+00_R;
        pts[ 1] = -9.4489927222288222340758013830322e-01_R;
        pts[ 2] = -8.1927932164400667834864158171690e-01_R;
        pts[ 3] = -6.3287615303186067766240485444366e-01_R;
        pts[ 4] = -3.9953094096534893226434979156697e-01_R;
        pts[ 5] = -1.3655293285492755486406185573969e-01_R;
        pts[ 6] = -pts[5];
        pts[ 7] = -pts[4];
        pts[ 8] = -pts[3];
        pts[ 9] = -pts[2];
        pts[10] = -pts[1];
        pts[11] = -pts[0];

        wts[ 0] = 1.5151515151515151515151515151515e-02_R;
        wts[ 1] = 9.1684517413196130668342594134079e-02_R;
        wts[ 2] = 1.5797470556437011516467106270034e-01_R;
        wts[ 3] = 2.1250841776102114535830207736687e-01_R;
        wts[ 4] = 2.5127560319920128029324441214760e-01_R;
        wts[ 5] = 2.7140524091069617700028833849960e-01_R;
        wts[ 6] = wts[5];
        wts[ 7] = wts[4];
        wts[ 8] = wts[3];
        wts[ 9] = wts[2];
        wts[10] = wts[1];
        wts[11] = wts[0];

        return;
      }

      // 2*13-3 = 23
    case TWENTYSECOND:
    case TWENTYTHIRD:
      {
        pts.resize (13);
        wts.resize(13);

        pts[ 0] = -1.0000000000000000000000000000000e+00_R;
        pts[ 1] = -9.5330984664216391189690546475545e-01_R;
        pts[ 2] = -8.4634756465187231686592560709875e-01_R;
        pts[ 3] = -6.8618846908175742607275903956636e-01_R;
        pts[ 4] = -4.8290982109133620174693723363693e-01_R;
        pts[ 5] = -2.4928693010623999256867370037423e-01_R;
        pts[ 6] = 0.;
        pts[ 7] = -pts[5];
        pts[ 8] = -pts[4];
        pts[ 9] = -pts[3];
        pts[10] = -pts[2];
        pts[11] = -pts[1];
        pts[12] = -pts[0];

        wts[ 0] = 1.2820512820512820512820512820513e-02_R;
        wts[ 1] = 7.7801686746818927793588988333134e-02_R;
        wts[ 2] = 1.3498192668960834911991476258937e-01_R;
        wts[ 3] = 1.8364686520355009200749425874681e-01_R;
        wts[ 4] = 2.2076779356611008608553400837940e-01_R;
        wts[ 5] = 2.4401579030667635645857814836016e-01_R;
        wts[ 6] = 2.5193084933344673604413864154124e-01_R;
        wts[ 7] = wts[5];
        wts[ 8] = wts[4];
        wts[ 9] = wts[3];
        wts[10] = wts[2];
        wts[11] = wts[1];
        wts[12] = wts[0];

        return;
      }

      // 2*14-3 = 25
    case TWENTYFOURTH:
    case TWENTYFIFTH:
      {
        pts.resize (14);
        wts.resize(14);

        pts[ 0] = -1.0000000000000000000000000000000e+00_R;
        pts[ 1] = -9.5993504526726090135510016201542e-01_R;
        pts[ 2] = -8.6780105383034725100022020290826e-01_R;
        pts[ 3] = -7.2886859909132614058467240052088e-01_R;
        pts[ 4] = -5.5063940292864705531662270585908e-01_R;
        pts[ 5] = -3.4272401334271284504390340364167e-01_R;
        pts[ 6] = -1.1633186888370386765877670973616e-01_R;
        pts[ 7] = -pts[6];
        pts[ 8] = -pts[5];
        pts[ 9] = -pts[4];
        pts[10] = -pts[3];
        pts[11] = -pts[2];
        pts[12] = -pts[1];
        pts[13] = -pts[0];

        wts[ 0] = 1.0989010989010989010989010989011e-02_R;
        wts[ 1] = 6.6837284497681284634070660746053e-02_R;
        wts[ 2] = 1.1658665589871165154099667065465e-01_R;
        wts[ 3] = 1.6002185176295214241282099798759e-01_R;
        wts[ 4] = 1.9482614937341611864033177837588e-01_R;
        wts[ 5] = 2.1912625300977075487116252395417e-01_R;
        wts[ 6] = 2.3161279446845705888962835729264e-01_R;
        wts[ 7] = wts[6];
        wts[ 8] = wts[5];
        wts[ 9] = wts[4];
        wts[10] = wts[3];
        wts[11] = wts[2];
        wts[12] = wts[1];
        wts[13] = wts[0];

        return;
      }

      // 2*15-3 = 27
    case TWENTYSIXTH:
    case TWENTYSEVENTH:
      {
        pts.resize (15);
        wts.resize(15);

        pts[ 0] = -1.0000000000000000000000000000000e+00_R;
        pts[ 1] = -9.6524592650383857279585139206960e-01_R;
        pts[ 2] = -8.8508204422297629882540163148223e-01_R;
        pts[ 3] = -7.6351968995181520070411847597629e-01_R;
        pts[ 4] = -6.0625320546984571112352993863673e-01_R;
        pts[ 5] = -4.2063805471367248092189693873858e-01_R;
        pts[ 6] = -2.1535395536379423822567944627292e-01_R;
        pts[ 7] = 0.;
        pts[ 8] = -pts[6];
        pts[ 9] = -pts[5];
        pts[10] = -pts[4];
        pts[11] = -pts[3];
        pts[12] = -pts[2];
        pts[13] = -pts[1];
        pts[14] = -pts[0];

        wts[ 0] = 9.5238095238095238095238095238095e-03_R;
        wts[ 1] = 5.8029893028601249096880584025282e-02_R;
        wts[ 2] = 1.0166007032571806760366617078880e-01_R;
        wts[ 3] = 1.4051169980242810946044680564367e-01_R;
        wts[ 4] = 1.7278964725360094905207709940835e-01_R;
        wts[ 5] = 1.9698723596461335609250034650741e-01_R;
        wts[ 6] = 2.1197358592682092012743007697722e-01_R;
        wts[ 7] = 2.1704811634881564951495021425091e-01_R;
        wts[ 8] = wts[6];
        wts[ 9] = wts[5];
        wts[10] = wts[4];
        wts[11] = wts[3];
        wts[12] = wts[2];
        wts[13] = wts[1];
        wts[14] = wts[0];

        return;
      }

      // 2*16-3 = 29
    case TWENTYEIGHTH:
    case TWENTYNINTH:
      {
        pts.resize (16);
        wts.resize(16);

        pts[ 0] = -1.0000000000000000000000000000000e+00_R;
        pts[ 1] = -9.6956804627021793295224273836746e-01_R;
        pts[ 2] = -8.9920053309347209299462826151985e-01_R;
        pts[ 3] = -7.9200829186181506393108827096315e-01_R;
        pts[ 4] = -6.5238870288249308946788321964058e-01_R;
        pts[ 5] = -4.8605942188713761178189078584687e-01_R;
        pts[ 6] = -2.9983046890076320809835345472230e-01_R;
        pts[ 7] = -1.0132627352194944784303300504592e-01_R;
        pts[ 8] = -pts[7];
        pts[ 9] = -pts[6];
        pts[10] = -pts[5];
        pts[11] = -pts[4];
        pts[12] = -pts[3];
        pts[13] = -pts[2];
        pts[14] = -pts[1];
        pts[15] = -pts[0];

        wts[ 0] = 8.3333333333333333333333333333333e-03_R;
        wts[ 1] = 5.0850361005919905403244919565455e-02_R;
        wts[ 2] = 8.9393697325930800991052080166084e-02_R;
        wts[ 3] = 1.2425538213251409834953633265731e-01_R;
        wts[ 4] = 1.5402698080716428081564494048499e-01_R;
        wts[ 5] = 1.7749191339170412530107566952836e-01_R;
        wts[ 6] = 1.9369002382520358431691359885352e-01_R;
        wts[ 7] = 2.0195830817822987148919912541094e-01_R;
        wts[ 8] = wts[7];
        wts[ 9] = wts[6];
        wts[10] = wts[5];
        wts[11] = wts[4];
        wts[12] = wts[3];
        wts[13] = wts[2];
        wts[14] = wts[1];
        wts[15] = wts[0];

        return;
      }

      // 2*17-3 = 31
    case THIRTIETH:
    case THIRTYFIRST:
      {
        pts.resize (17);
        wts.resize(17);

        pts[ 0] = -1.0000000000000000000000000000000e+00_R;
        pts[ 1] = -9.7313217663141831415697950187372e-01_R;
        pts[ 2] = -9.1087999591557359562380250639773e-01_R;
        pts[ 3] = -8.1569625122177030710675055323753e-01_R;
        pts[ 4] = -6.9102898062768470539491935737245e-01_R;
        pts[ 5] = -5.4138539933010153912373340750406e-01_R;
        pts[ 6] = -3.7217443356547704190723468073526e-01_R;
        pts[ 7] = -1.8951197351831738830426301475311e-01_R;
        pts[ 8] = 0.;
        pts[ 9] = -pts[7];
        pts[10] = -pts[6];
        pts[11] = -pts[5];
        pts[12] = -pts[4];
        pts[13] = -pts[3];
        pts[14] = -pts[2];
        pts[15] = -pts[1];
        pts[16] = -pts[0];

        wts[ 0] = 7.3529411764705882352941176470588e-03_R;
        wts[ 1] = 4.4921940543254209647400954623212e-02_R;
        wts[ 2] = 7.9198270503687119190264429952835e-02_R;
        wts[ 3] = 1.1059290900702816137577270522008e-01_R;
        wts[ 4] = 1.3798774620192655905620157495403e-01_R;
        wts[ 5] = 1.6039466199762153951632836586475e-01_R;
        wts[ 6] = 1.7700425351565787043694574536329e-01_R;
        wts[ 7] = 1.8721633967761923589208848286062e-01_R;
        wts[ 8] = 1.9066187475346943329940724702825e-01_R;
        wts[ 9] = wts[7];
        wts[10] = wts[6];
        wts[11] = wts[5];
        wts[12] = wts[4];
        wts[13] = wts[3];
        wts[14] = wts[2];
        wts[15] = wts[1];
        wts[16] = wts[0];

        return;
      }

      // 2*18-3 = 33
    case THIRTYSECOND:
    case THIRTYTHIRD:
      {
        pts.resize (18);
        wts.resize(18);

        pts[ 0] = -1.0000000000000000000000000000000e+00_R;
        pts[ 1] = -9.7610555741219854286451892434170e-01_R;
        pts[ 2] = -9.2064918534753387383785462543128e-01_R;
        pts[ 3] = -8.3559353521809021371364636232794e-01_R;
        pts[ 4] = -7.2367932928324268130621036530207e-01_R;
        pts[ 5] = -5.8850483431866176117353589319356e-01_R;
        pts[ 6] = -4.3441503691212397534228713674067e-01_R;
        pts[ 7] = -2.6636265287828098416766533202560e-01_R;
        pts[ 8] = -8.9749093484652111022645010088562e-02_R;
        pts[ 9] = -pts[8];
        pts[10] = -pts[7];
        pts[11] = -pts[6];
        pts[12] = -pts[5];
        pts[13] = -pts[4];
        pts[14] = -pts[3];
        pts[15] = -pts[2];
        pts[16] = -pts[1];
        pts[17] = -pts[0];

        wts[ 0] = 6.5359477124183006535947712418301e-03_R;
        wts[ 1] = 3.9970628810914066137599176410101e-02_R;
        wts[ 2] = 7.0637166885633664999222960167786e-02_R;
        wts[ 3] = 9.9016271717502802394423605318672e-02_R;
        wts[ 4] = 1.2421053313296710026339635889675e-01_R;
        wts[ 5] = 1.4541196157380226798300321049443e-01_R;
        wts[ 6] = 1.6193951723760248926432670670023e-01_R;
        wts[ 7] = 1.7326210948945622601061440382668e-01_R;
        wts[ 8] = 1.7901586343970308229381880694353e-01_R;
        wts[ 9] = wts[8];
        wts[10] = wts[7];
        wts[11] = wts[6];
        wts[12] = wts[5];
        wts[13] = wts[4];
        wts[14] = wts[3];
        wts[15] = wts[2];
        wts[16] = wts[1];
        wts[17] = wts[0];

        return;
      }

      // 2*19-3 = 35
    case THIRTYFOURTH:
    case THIRTYFIFTH:
      {
        pts.resize (19);
        wts.resize(19);

        pts[ 0] = -1.0000000000000000000000000000000e+00_R;
        pts[ 1] = -9.7861176622208009515263406311022e-01_R;
        pts[ 2] = -9.2890152815258624371794025879655e-01_R;
        pts[ 3] = -8.5246057779664609308595597004106e-01_R;
        pts[ 4] = -7.5149420255261301416363748963394e-01_R;
        pts[ 5] = -6.2890813726522049776683230622873e-01_R;
        pts[ 6] = -4.8822928568071350277790963762492e-01_R;
        pts[ 7] = -3.3350484782449861029850010384493e-01_R;
        pts[ 8] = -1.6918602340928157137515415344488e-01_R;
        pts[ 9] = 0.;
        pts[10] = -pts[8];
        pts[11] = -pts[7];
        pts[12] = -pts[6];
        pts[13] = -pts[5];
        pts[14] = -pts[4];
        pts[15] = -pts[3];
        pts[16] = -pts[2];
        pts[17] = -pts[1];
        pts[18] = -pts[0];

        wts[ 0] = 5.8479532163742690058479532163743e-03_R;
        wts[ 1] = 3.5793365186176477115425569035122e-02_R;
        wts[ 2] = 6.3381891762629736851695690418317e-02_R;
        wts[ 3] = 8.9131757099207084448008790556153e-02_R;
        wts[ 4] = 1.1231534147730504407091001546378e-01_R;
        wts[ 5] = 1.3226728044875077692604673390973e-01_R;
        wts[ 6] = 1.4841394259593888500968064366841e-01_R;
        wts[ 7] = 1.6029092404406124197991096818359e-01_R;
        wts[ 8] = 1.6755658452714286727013727774026e-01_R;
        wts[ 9] = 1.7000191928482723464467271561652e-01_R;
        wts[10] = wts[8];
        wts[11] = wts[7];
        wts[12] = wts[6];
        wts[13] = wts[5];
        wts[14] = wts[4];
        wts[15] = wts[3];
        wts[16] = wts[2];
        wts[17] = wts[1];
        wts[18] = wts[0];

        return;
      }

      // 2*20-3 = 37
    case THIRTYSIXTH:
    case THIRTYSEVENTH:
      {
        pts.resize (20);
        wts.resize(20);

        pts[ 0] = -1.0000000000000000000000000000000e+00_R;
        pts[ 1] = -9.8074370489391417192544643858423e-01_R;
        pts[ 2] = -9.3593449881266543571618158493063e-01_R;
        pts[ 3] = -8.6687797808995014130984721461629e-01_R;
        pts[ 4] = -7.7536826095205587041431752759469e-01_R;
        pts[ 5] = -6.6377640229031128984640332297116e-01_R;
        pts[ 6] = -5.3499286403188626164813596182898e-01_R;
        pts[ 7] = -3.9235318371390929938647470381582e-01_R;
        pts[ 8] = -2.3955170592298649518240135692709e-01_R;
        pts[ 9] = -8.0545937238821837975944518159554e-02_R;
        pts[10] = -pts[9];
        pts[11] = -pts[8];
        pts[12] = -pts[7];
        pts[13] = -pts[6];
        pts[14] = -pts[5];
        pts[15] = -pts[4];
        pts[16] = -pts[3];
        pts[17] = -pts[2];
        pts[18] = -pts[1];
        pts[19] = -pts[0];

        wts[ 0] = 5.2631578947368421052631578947368e-03_R;
        wts[ 1] = 3.2237123188488941491605028117294e-02_R;
        wts[ 2] = 5.7181802127566826004753627173243e-02_R;
        wts[ 3] = 8.0631763996119603144776846113721e-02_R;
        wts[ 4] = 1.0199149969945081568378120573289e-01_R;
        wts[ 5] = 1.2070922762867472509942970500239e-01_R;
        wts[ 6] = 1.3630048235872418448978079298903e-01_R;
        wts[ 7] = 1.4836155407091682581471301373397e-01_R;
        wts[ 8] = 1.5658010264747548715816989679364e-01_R;
        wts[ 9] = 1.6074328638784574900772672644908e-01_R;
        wts[10] = wts[9];
        wts[11] = wts[8];
        wts[12] = wts[7];
        wts[13] = wts[6];
        wts[14] = wts[5];
        wts[15] = wts[4];
        wts[16] = wts[3];
        wts[17] = wts[2];
        wts[18] = wts[1];
        wts[19] = wts[0];

        return;
      }

      // 2*21-3 = 39
    case THIRTYEIGHTH:
    case THIRTYNINTH:
      {
        pts.resize (21);
        wts.resize(21);

        pts[ 0] = -1.0000000000000000000000000000000e+00_R;
        pts[ 1] = -9.8257229660454802823448127655541e-01_R;
        pts[ 2] = -9.4197629695974553429610265066144e-01_R;
        pts[ 3] = -8.7929475532359046445115359630494e-01_R;
        pts[ 4] = -7.9600192607771240474431258966036e-01_R;
        pts[ 5] = -6.9405102606222323262731639319467e-01_R;
        pts[ 6] = -5.7583196026183068692702187033809e-01_R;
        pts[ 7] = -4.4411578327900210119451634960735e-01_R;
        pts[ 8] = -3.0198985650876488727535186785875e-01_R;
        pts[ 9] = -1.5278551580218546600635832848567e-01_R;
        pts[10] = 0.;
        pts[11] = -pts[9];
        pts[12] = -pts[8];
        pts[13] = -pts[7];
        pts[14] = -pts[6];
        pts[15] = -pts[5];
        pts[16] = -pts[4];
        pts[17] = -pts[3];
        pts[18] = -pts[2];
        pts[19] = -pts[1];
        pts[20] = -pts[0];

        wts[ 0] = 4.7619047619047619047619047619048e-03_R;
        wts[ 1] = 2.9184840098505458609458543613171e-02_R;
        wts[ 2] = 5.1843169000849625072722971852830e-02_R;
        wts[ 3] = 7.3273918185074144252547861041894e-02_R;
        wts[ 4] = 9.2985467957886065301137664149214e-02_R;
        wts[ 5] = 1.1051708321912333526700048678439e-01_R;
        wts[ 6] = 1.2545812119086894801515753570800e-01_R;
        wts[ 7] = 1.3745846286004134358089961741515e-01_R;
        wts[ 8] = 1.4623686244797745926727053063439e-01_R;
        wts[ 9] = 1.5158757511168138445325068150529e-01_R;
        wts[10] = 1.5338519033217494855158440506754e-01_R;
        wts[11] = wts[9];
        wts[12] = wts[8];
        wts[13] = wts[7];
        wts[14] = wts[6];
        wts[15] = wts[5];
        wts[16] = wts[4];
        wts[17] = wts[3];
        wts[18] = wts[2];
        wts[19] = wts[1];
        wts[20] = wts[0];

        return;
      }

      // 2*22-3 = 41
    case FORTIETH:
    case FORTYFIRST:
      {
        pts.resize (22);
        wts.resize(22);

        pts[ 0] = -1.0000000000000000000000000000000e+00_R;
        pts[ 1] = -9.8415243845764617655228962221207e-01_R;
        pts[ 2] = -9.4720428399922868052421376661573e-01_R;
        pts[ 3] = -8.9006229019090447052965782577909e-01_R;
        pts[ 4] = -8.1394892761192113604544184805614e-01_R;
        pts[ 5] = -7.2048723996120215811988189639847e-01_R;
        pts[ 6] = -6.1166943828425897122621160586993e-01_R;
        pts[ 7] = -4.8981487518990234980875123568327e-01_R;
        pts[ 8] = -3.5752071013891953806095728024018e-01_R;
        pts[ 9] = -2.1760658515928504178795509346539e-01_R;
        pts[10] = -7.3054540010898334761088790464107e-02_R;
        pts[11] = -pts[10];
        pts[12] = -pts[9];
        pts[13] = -pts[8];
        pts[14] = -pts[7];
        pts[15] = -pts[6];
        pts[16] = -pts[5];
        pts[17] = -pts[4];
        pts[18] = -pts[3];
        pts[19] = -pts[2];
        pts[20] = -pts[1];
        pts[21] = -pts[0];

        wts[ 0] = 4.3290043290043290043290043290043e-03_R;
        wts[ 1] = 2.6545747682501757911627904520543e-02_R;
        wts[ 2] = 4.7214465293740752123775734864792e-02_R;
        wts[ 3] = 6.6865605864553076012404194157097e-02_R;
        wts[ 4] = 8.5090060391838447815711236095748e-02_R;
        wts[ 5] = 1.0150057480164767437243730374960e-01_R;
        wts[ 6] = 1.1574764465393906659003636772146e-01_R;
        wts[ 7] = 1.2752769665343027553084445930883e-01_R;
        wts[ 8] = 1.3658968861374142668617736220617e-01_R;
        wts[ 9] = 1.4274049227136140033623599356679e-01_R;
        wts[10] = 1.4584901944424179361642043947997e-01_R;
        wts[11] = wts[10];
        wts[12] = wts[9];
        wts[13] = wts[8];
        wts[14] = wts[7];
        wts[15] = wts[6];
        wts[16] = wts[5];
        wts[17] = wts[4];
        wts[18] = wts[3];
        wts[19] = wts[2];
        wts[20] = wts[1];
        wts[21] = wts[0];

        return;
      }

      // 2*23-3 = 43
    case FORTYSECOND:
    case FORTYTHIRD:
      {
        pts.resize (23);
        wts.resize(23);

        pts[ 0] = -1.0000000000000000000000000000000e+00_R;
        pts[ 1] = -9.8552715587873257808146276673810e-01_R;
        pts[ 2] = -9.5175795571071020413563967985143e-01_R;
        pts[ 3] = -8.9945855804034501095016032034737e-01_R;
        pts[ 4] = -8.2965109665128588622320061929000e-01_R;
        pts[ 5] = -7.4369504117206068394516354306700e-01_R;
        pts[ 6] = -6.4326364446013620847614553360277e-01_R;
        pts[ 7] = -5.3031177113684416813011532015230e-01_R;
        pts[ 8] = -4.0703793791447482919595048821510e-01_R;
        pts[ 9] = -2.7584154894579306710687763267914e-01_R;
        pts[10] = -1.3927620404066839859186261298277e-01_R;
        pts[11] = 0.;
        pts[12] = -pts[10];
        pts[13] = -pts[9];
        pts[14] = -pts[8];
        pts[15] = -pts[7];
        pts[16] = -pts[6];
        pts[17] = -pts[5];
        pts[18] = -pts[4];
        pts[19] = -pts[3];
        pts[20] = -pts[2];
        pts[21] = -pts[1];
        pts[22] = -pts[0];

        wts[ 0] = 3.9525691699604743083003952569170e-03_R;
        wts[ 1] = 2.4248600771531736517399658937097e-02_R;
        wts[ 2] = 4.3175871170241834748876465612042e-02_R;
        wts[ 3] = 6.1252477129554206381382847440355e-02_R;
        wts[ 4] = 7.8135449475569989741934255347965e-02_R;
        wts[ 5] = 9.3497246163512341833500706906697e-02_R;
        wts[ 6] = 1.0703910172433651153518362791547e-01_R;
        wts[ 7] = 1.1849751066274913130212600472426e-01_R;
        wts[ 8] = 1.2764947470175887663614855305567e-01_R;
        wts[ 9] = 1.3431687263860381990156489770071e-01_R;
        wts[10] = 1.3836993638580739452350273386294e-01_R;
        wts[11] = 1.3972978001274736514015970647975e-01_R;
        wts[12] = wts[10];
        wts[13] = wts[9];
        wts[14] = wts[8];
        wts[15] = wts[7];
        wts[16] = wts[6];
        wts[17] = wts[5];
        wts[18] = wts[4];
        wts[19] = wts[3];
        wts[20] = wts[2];
        wts[21] = wts[1];
        wts[22] = wts[0];

        return;
      }

    default:
      libmesh_error_msg("Quadrature rule " << order << " not supported!");
    }
}

/// Every 1D rule the table above holds, built once and shared by all callers
struct Table1D
{
  std::array<std::vector<Real>, QGaussLobatto::max_points_1D + 1> points, weights;

  Table1D()
  {
    for (unsigned int n = 2; n <= QGaussLobatto::max_points_1D; ++n)
      fill_1D(static_cast<Order>(2 * n - 3), points[n], weights[n]);
  }
};

const Table1D & table_1D()
{
  static const Table1D table;
  return table;
}

} // anonymous namespace


unsigned int QGaussLobatto::n_points_1D(const Order order)
{
  // The table pairs the orders 2n-4 and 2n-3 with the same n points, and there is no
  // one-point rule, so the lowest orders share the two-point rule.
  return static_cast<unsigned int>(order) / 2 + 2;
}


const std::vector<Real> & QGaussLobatto::points_1D(const unsigned int n)
{
  libmesh_error_msg_if(n < 2 || n > max_points_1D,
                       "No " << n << " point Gauss-Lobatto rule is tabulated; "
                       "the rules run from 2 to " << max_points_1D << " points.");

  return table_1D().points[n];
}


const std::vector<Real> & QGaussLobatto::weights_1D(const unsigned int n)
{
  libmesh_error_msg_if(n < 2 || n > max_points_1D,
                       "No " << n << " point Gauss-Lobatto rule is tabulated; "
                       "the rules run from 2 to " << max_points_1D << " points.");

  return table_1D().weights[n];
}


void QGaussLobatto::init_1D()
{
  const unsigned int n = n_points_1D(get_order());

  libmesh_error_msg_if(n > max_points_1D,
                       "Quadrature rule " << get_order() << " not supported!");

  const std::vector<Real> & pts = points_1D(n);
  const std::vector<Real> & wts = weights_1D(n);

  _points.resize (n);
  _weights.resize(n);

  for (const auto i : index_range(pts))
    {
      _points[i](0) = pts[i];
      _weights[i]   = wts[i];
    }
}

} // namespace libMesh
