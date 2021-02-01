/** @file
 *****************************************************************************
 * @author     This file is part of libff, developed by SCIPR Lab
 *             and contributors (see AUTHORS).
 * @copyright  MIT license (see LICENSE file)
 *****************************************************************************/

#include <libff/algebra/curves/edwards/edwards_fields.hpp>

namespace libff {

bigint<edwards_r_limbs> edwards_modulus_r;
bigint<edwards_q_limbs> edwards_modulus_q;

void init_edwards_fields()
{
    typedef bigint<edwards_r_limbs> bigint_r;
    typedef bigint<edwards_q_limbs> bigint_q;

    assert(sizeof(mp_limb_t) == 8 || sizeof(mp_limb_t) == 4); // Montgomery assumes this

    /* parameters for scalar field Fr */

    edwards_modulus_r = bigint_r("1552511030102430251236801561344621993261920897571225601");
    assert(edwards_Fr::modulus_is_valid());
    if (sizeof(mp_limb_t) == 8)
    {
        edwards_Fr::Rsquared = bigint_r("621738487827897760168419760282818735947979812540885779");
        edwards_Fr::Rcubed = bigint_r("899968968216802386013510389846941393831065658679774050");
        edwards_Fr::inv = 0xdde553277fffffff;
    }
    if (sizeof(mp_limb_t) == 4)
    {
        edwards_Fr::Rsquared = bigint_r("621738487827897760168419760282818735947979812540885779");
        edwards_Fr::Rcubed = bigint_r("899968968216802386013510389846941393831065658679774050");
        edwards_Fr::inv = 0x7fffffff;
    }
    edwards_Fr::num_bits = 181;
    edwards_Fr::euler = bigint_r("776255515051215125618400780672310996630960448785612800");
    edwards_Fr::s = 31;
    edwards_Fr::t = bigint_r("722944284836962004768104088187507350585386575");
    edwards_Fr::t_minus_1_over_2 = bigint_r("361472142418481002384052044093753675292693287");
    edwards_Fr::multiplicative_generator = edwards_Fr("19");
    edwards_Fr::root_of_unity = edwards_Fr("695314865466598274460565335217615316274564719601897184");
    edwards_Fr::nqr = edwards_Fr("11");
    edwards_Fr::nqr_to_t = edwards_Fr("1326707053668679463752768729767248251415639579872144553");

    /* parameters for base field Fq */

    edwards_modulus_q = bigint_q("6210044120409721004947206240885978274523751269793792001");
    assert(edwards_Fq::modulus_is_valid());
    if (sizeof(mp_limb_t) == 8)
    {
        edwards_Fq::Rsquared = bigint_q("5943559676554581037560514598978484097352477055348195432");
        edwards_Fq::Rcubed = bigint_q("1081560488703514202058739223469726982199727506489234349");
        edwards_Fq::inv = 0x76eb690b7fffffff;
    }
    if (sizeof(mp_limb_t) == 4)
    {
        edwards_Fq::Rsquared = bigint_q("5943559676554581037560514598978484097352477055348195432");
        edwards_Fq::Rcubed = bigint_q("1081560488703514202058739223469726982199727506489234349");
        edwards_Fq::inv = 0x7fffffff;
    }
    edwards_Fq::num_bits = 183;
    edwards_Fq::euler = bigint_q("3105022060204860502473603120442989137261875634896896000");
    edwards_Fq::s = 31;
    edwards_Fq::t = bigint_q("2891777139347848019072416350658041552884388375");
    edwards_Fq::t_minus_1_over_2 = bigint_q("1445888569673924009536208175329020776442194187");
    edwards_Fq::multiplicative_generator = edwards_Fq("61");
    edwards_Fq::root_of_unity = edwards_Fq("4692813029219384139894873043933463717810008194158530536");
    edwards_Fq::nqr = edwards_Fq("23");
    edwards_Fq::nqr_to_t = edwards_Fq("2626736066325740702418554487368721595489070118548299138");

    /* parameters for twist field Fq3 */

    edwards_Fq3::euler = bigint<3*edwards_q_limbs>("119744082713971502962992613191067836698205043373978948903839934564152994858051284658545502971203325031831647424413111161318314144765646525057914792711854057586688000");
    edwards_Fq3::s = 31;
    edwards_Fq3::t = bigint<3*edwards_q_limbs>("111520367408144756185815309352304634357062208814526860512643991563611659089151103662834971185031649686239331424621037357783237607000066456438894190557165125");
    edwards_Fq3::t_minus_1_over_2 = bigint<3*edwards_q_limbs>("55760183704072378092907654676152317178531104407263430256321995781805829544575551831417485592515824843119665712310518678891618803500033228219447095278582562");
    edwards_Fq3::non_residue = edwards_Fq("61");
    edwards_Fq3::nqr = edwards_Fq3(edwards_Fq("23"),edwards_Fq("0"),edwards_Fq("0"));
    edwards_Fq3::nqr_to_t = edwards_Fq3(edwards_Fq("104810943629412208121981114244673004633270996333237516"),edwards_Fq("0"),edwards_Fq("0"));
    edwards_Fq3::Frobenius_coeffs_c1[0] = edwards_Fq("1");
    edwards_Fq3::Frobenius_coeffs_c1[1] = edwards_Fq("1073752683758513276629212192812154536507607213288832061");
    edwards_Fq3::Frobenius_coeffs_c1[2] = edwards_Fq("5136291436651207728317994048073823738016144056504959939");
    edwards_Fq3::Frobenius_coeffs_c2[0] = edwards_Fq("1");
    edwards_Fq3::Frobenius_coeffs_c2[1] = edwards_Fq("5136291436651207728317994048073823738016144056504959939");
    edwards_Fq3::Frobenius_coeffs_c2[2] = edwards_Fq("1073752683758513276629212192812154536507607213288832061");

    /* parameters for Fq6 */

    edwards_Fq6::euler = bigint<6*edwards_q_limbs>("286772906900208978053659358438511744950463073927050691441013733733308852103113703187330276692859302954336930908810347208181904968549176866575351093481768308610598000"
                                                   "12474722471292022005286908893541416599742789523143737959708214461870605874121090272322517073698670457333498341250241297120220551099521678676677016084977733861376000");
    edwards_Fq6::s = 32;
    edwards_Fq6::t = bigint<6*edwards_q_limbs>("133539040992133774819625242817453921283810899558966370039160239386050427380536544689596759478030401726410226388280441615815230164932635084993835116495270425738113932"
                                               "36962483483968312854546545489477653335660132357451942730559136002264267706368626222164610900775841650242683636606783278808100405661166495958658408413165125");
    edwards_Fq6::t_minus_1_over_2 = bigint<6*edwards_q_limbs>("667695204960668874098126214087269606419054497794831850195801196930252136902682723447983797390152008632051131941402208079076150824663175424969175582476352128690569661"
                                                              "8481241741984156427273272744738826667830066178725971365279568001132133853184313111082305450387920825121341818303391639404050202830583247979329204206582562");
    edwards_Fq6::non_residue = edwards_Fq("61");
    edwards_Fq6::nqr = edwards_Fq6(edwards_Fq3(edwards_Fq("5"),edwards_Fq("0"),edwards_Fq("0")),edwards_Fq3::one());
    edwards_Fq6::nqr_to_t = edwards_Fq6(edwards_Fq3::zero(),edwards_Fq3(edwards_Fq("0"),edwards_Fq("6018622460271751604575462891699668290753365582464183006"),edwards_Fq("0")));
    edwards_Fq6::Frobenius_coeffs_c1[0] = edwards_Fq("1");
    edwards_Fq6::Frobenius_coeffs_c1[1] = edwards_Fq("1073752683758513276629212192812154536507607213288832062");
    edwards_Fq6::Frobenius_coeffs_c1[2] = edwards_Fq("1073752683758513276629212192812154536507607213288832061");
    edwards_Fq6::Frobenius_coeffs_c1[3] = edwards_Fq("6210044120409721004947206240885978274523751269793792000");
    edwards_Fq6::Frobenius_coeffs_c1[4] = edwards_Fq("5136291436651207728317994048073823738016144056504959939");
    edwards_Fq6::Frobenius_coeffs_c1[5] = edwards_Fq("5136291436651207728317994048073823738016144056504959940");
    edwards_Fq6::my_Fp2::non_residue = edwards_Fq3::non_residue;
}

} // namespace libff
