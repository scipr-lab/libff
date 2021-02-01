/** @file
 *****************************************************************************
 * @author     This file is part of libff, developed by SCIPR Lab
 *             and contributors (see AUTHORS).
 * @copyright  MIT license (see LICENSE file)
 *****************************************************************************/

#include <libff/algebra/curves/bn128/bn128_fields.hpp>

namespace libff {

bigint<bn128_r_limbs> bn128_modulus_r;
bigint<bn128_q_limbs> bn128_modulus_q;

void init_bn128_fields()
{
    bn::Param::init(); // init ate-pairing library

    typedef bigint<bn128_r_limbs> bigint_r;
    typedef bigint<bn128_q_limbs> bigint_q;

    assert(sizeof(mp_limb_t) == 8 || sizeof(mp_limb_t) == 4); // Montgomery assumes this

    /* parameters for scalar field Fr */
    bn128_modulus_r = bigint_r("21888242871839275222246405745257275088548364400416034343698204186575808495617");
    assert(bn128_Fr::modulus_is_valid());
    if (sizeof(mp_limb_t) == 8)
    {
        bn128_Fr::Rsquared = bigint_r("944936681149208446651664254269745548490766851729442924617792859073125903783");
        bn128_Fr::Rcubed = bigint_r("5866548545943845227489894872040244720403868105578784105281690076696998248512");
        bn128_Fr::inv = 0xc2e1f593efffffff;
    }
    if (sizeof(mp_limb_t) == 4)
    {
        bn128_Fr::Rsquared = bigint_r("944936681149208446651664254269745548490766851729442924617792859073125903783");
        bn128_Fr::Rcubed = bigint_r("5866548545943845227489894872040244720403868105578784105281690076696998248512");
        bn128_Fr::inv = 0xefffffff;
    }
    bn128_Fr::num_bits = 254;
    bn128_Fr::euler = bigint_r("10944121435919637611123202872628637544274182200208017171849102093287904247808");
    bn128_Fr::s = 28;
    bn128_Fr::t = bigint_r("81540058820840996586704275553141814055101440848469862132140264610111");
    bn128_Fr::t_minus_1_over_2 = bigint_r("40770029410420498293352137776570907027550720424234931066070132305055");
    bn128_Fr::multiplicative_generator = bn128_Fr("5");
    bn128_Fr::root_of_unity = bn128_Fr("19103219067921713944291392827692070036145651957329286315305642004821462161904");
    bn128_Fr::nqr = bn128_Fr("5");
    bn128_Fr::nqr_to_t = bn128_Fr("19103219067921713944291392827692070036145651957329286315305642004821462161904");

    /* parameters for base field Fq */
    bn128_modulus_q = bigint_q("21888242871839275222246405745257275088696311157297823662689037894645226208583");
    assert(bn128_Fq::modulus_is_valid());
    if (sizeof(mp_limb_t) == 8)
    {
        bn128_Fq::Rsquared = bigint_q("3096616502983703923843567936837374451735540968419076528771170197431451843209");
        bn128_Fq::Rcubed = bigint_q("14921786541159648185948152738563080959093619838510245177710943249661917737183");
        bn128_Fq::inv = 0x87d20782e4866389;
    }
    if (sizeof(mp_limb_t) == 4)
    {
        bn128_Fq::Rsquared = bigint_q("3096616502983703923843567936837374451735540968419076528771170197431451843209");
        bn128_Fq::Rcubed = bigint_q("14921786541159648185948152738563080959093619838510245177710943249661917737183");
        bn128_Fq::inv = 0xe4866389;
    }
    bn128_Fq::num_bits = 254;
    bn128_Fq::euler = bigint_q("10944121435919637611123202872628637544348155578648911831344518947322613104291");
    bn128_Fq::s = 1;
    bn128_Fq::t = bigint_q("10944121435919637611123202872628637544348155578648911831344518947322613104291");
    bn128_Fq::t_minus_1_over_2 = bigint_q("5472060717959818805561601436314318772174077789324455915672259473661306552145");
    bn128_Fq::multiplicative_generator = bn128_Fq("3");
    bn128_Fq::root_of_unity = bn128_Fq("21888242871839275222246405745257275088696311157297823662689037894645226208582");
    bn128_Fq::nqr = bn128_Fq("3");
    bn128_Fq::nqr_to_t = bn128_Fq("21888242871839275222246405745257275088696311157297823662689037894645226208582");
}

} // namespace libff
