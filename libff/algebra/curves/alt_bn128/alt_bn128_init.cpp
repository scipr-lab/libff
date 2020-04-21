/** @file
 *****************************************************************************
 * @author     This file is part of libff, developed by SCIPR Lab
 *             and contributors (see AUTHORS).
 * @copyright  MIT license (see LICENSE file)
 *****************************************************************************/

#include <libff/algebra/curves/alt_bn128/alt_bn128_g1.hpp>
#include <libff/algebra/curves/alt_bn128/alt_bn128_g2.hpp>
#include <libff/algebra/curves/alt_bn128/alt_bn128_init.hpp>

namespace libff {

alt_bn128_Fq alt_bn128_coeff_b;
alt_bn128_Fq2 alt_bn128_twist;
alt_bn128_Fq2 alt_bn128_twist_coeff_b;
alt_bn128_Fq alt_bn128_twist_mul_by_b_c0;
alt_bn128_Fq alt_bn128_twist_mul_by_b_c1;
alt_bn128_Fq2 alt_bn128_twist_mul_by_q_X;
alt_bn128_Fq2 alt_bn128_twist_mul_by_q_Y;

bigint<alt_bn128_q_limbs> alt_bn128_ate_loop_count;
bool alt_bn128_ate_is_loop_count_neg;
bigint<12*alt_bn128_q_limbs> alt_bn128_final_exponent;
bigint<alt_bn128_q_limbs> alt_bn128_final_exponent_z;
bool alt_bn128_final_exponent_is_z_neg;

void init_alt_bn128_params()
{
    init_alt_bn128_fields();

    /* choice of short Weierstrass curve and its twist */

    alt_bn128_coeff_b = alt_bn128_Fq("3");
    alt_bn128_twist = alt_bn128_Fq2(alt_bn128_Fq("9"), alt_bn128_Fq("1"));
    alt_bn128_twist_coeff_b = alt_bn128_coeff_b * alt_bn128_twist.inverse();
    alt_bn128_twist_mul_by_b_c0 = alt_bn128_coeff_b * alt_bn128_Fq2::non_residue;
    alt_bn128_twist_mul_by_b_c1 = alt_bn128_coeff_b * alt_bn128_Fq2::non_residue;
    alt_bn128_twist_mul_by_q_X = alt_bn128_Fq2(alt_bn128_Fq("21575463638280843010398324269430826099269044274347216827212613867836435027261"),
                                           alt_bn128_Fq("10307601595873709700152284273816112264069230130616436755625194854815875713954"));
    alt_bn128_twist_mul_by_q_Y = alt_bn128_Fq2(alt_bn128_Fq("2821565182194536844548159561693502659359617185244120367078079554186484126554"),
                                           alt_bn128_Fq("3505843767911556378687030309984248845540243509899259641013678093033130930403"));

    /* choice of group G1 */
    alt_bn128_G1::G1_zero = alt_bn128_G1(alt_bn128_Fq::zero(),
                                     alt_bn128_Fq::one(),
                                     alt_bn128_Fq::zero());
    alt_bn128_G1::G1_one = alt_bn128_G1(alt_bn128_Fq("1"),
                                    alt_bn128_Fq("2"),
                                    alt_bn128_Fq::one());
    alt_bn128_G1::wnaf_window_table.resize(0);
    alt_bn128_G1::wnaf_window_table.push_back(11);
    alt_bn128_G1::wnaf_window_table.push_back(24);
    alt_bn128_G1::wnaf_window_table.push_back(60);
    alt_bn128_G1::wnaf_window_table.push_back(127);

    alt_bn128_G1::fixed_base_exp_window_table.resize(0);
    // window 1 is unbeaten in [-inf, 4.99]
    alt_bn128_G1::fixed_base_exp_window_table.push_back(1);
    // window 2 is unbeaten in [4.99, 10.99]
    alt_bn128_G1::fixed_base_exp_window_table.push_back(5);
    // window 3 is unbeaten in [10.99, 32.29]
    alt_bn128_G1::fixed_base_exp_window_table.push_back(11);
    // window 4 is unbeaten in [32.29, 55.23]
    alt_bn128_G1::fixed_base_exp_window_table.push_back(32);
    // window 5 is unbeaten in [55.23, 162.03]
    alt_bn128_G1::fixed_base_exp_window_table.push_back(55);
    // window 6 is unbeaten in [162.03, 360.15]
    alt_bn128_G1::fixed_base_exp_window_table.push_back(162);
    // window 7 is unbeaten in [360.15, 815.44]
    alt_bn128_G1::fixed_base_exp_window_table.push_back(360);
    // window 8 is unbeaten in [815.44, 2373.07]
    alt_bn128_G1::fixed_base_exp_window_table.push_back(815);
    // window 9 is unbeaten in [2373.07, 6977.75]
    alt_bn128_G1::fixed_base_exp_window_table.push_back(2373);
    // window 10 is unbeaten in [6977.75, 7122.23]
    alt_bn128_G1::fixed_base_exp_window_table.push_back(6978);
    // window 11 is unbeaten in [7122.23, 57818.46]
    alt_bn128_G1::fixed_base_exp_window_table.push_back(7122);
    // window 12 is never the best
    alt_bn128_G1::fixed_base_exp_window_table.push_back(0);
    // window 13 is unbeaten in [57818.46, 169679.14]
    alt_bn128_G1::fixed_base_exp_window_table.push_back(57818);
    // window 14 is never the best
    alt_bn128_G1::fixed_base_exp_window_table.push_back(0);
    // window 15 is unbeaten in [169679.14, 439758.91]
    alt_bn128_G1::fixed_base_exp_window_table.push_back(169679);
    // window 16 is unbeaten in [439758.91, 936073.41]
    alt_bn128_G1::fixed_base_exp_window_table.push_back(439759);
    // window 17 is unbeaten in [936073.41, 4666554.74]
    alt_bn128_G1::fixed_base_exp_window_table.push_back(936073);
    // window 18 is never the best
    alt_bn128_G1::fixed_base_exp_window_table.push_back(0);
    // window 19 is unbeaten in [4666554.74, 7580404.42]
    alt_bn128_G1::fixed_base_exp_window_table.push_back(4666555);
    // window 20 is unbeaten in [7580404.42, 34552892.20]
    alt_bn128_G1::fixed_base_exp_window_table.push_back(7580404);
    // window 21 is never the best
    alt_bn128_G1::fixed_base_exp_window_table.push_back(0);
    // window 22 is unbeaten in [34552892.20, inf]
    alt_bn128_G1::fixed_base_exp_window_table.push_back(34552892);

    /* choice of group G2 */

    alt_bn128_G2::G2_zero = alt_bn128_G2(alt_bn128_Fq2::zero(),
                                     alt_bn128_Fq2::one(),
                                     alt_bn128_Fq2::zero());

    alt_bn128_G2::G2_one = alt_bn128_G2(alt_bn128_Fq2(alt_bn128_Fq("10857046999023057135944570762232829481370756359578518086990519993285655852781"),
                                                alt_bn128_Fq("11559732032986387107991004021392285783925812861821192530917403151452391805634")),
                                    alt_bn128_Fq2(alt_bn128_Fq("8495653923123431417604973247489272438418190587263600148770280649306958101930"),
                                                alt_bn128_Fq("4082367875863433681332203403145435568316851327593401208105741076214120093531")),
                                    alt_bn128_Fq2::one());
    alt_bn128_G2::wnaf_window_table.resize(0);
    alt_bn128_G2::wnaf_window_table.push_back(5);
    alt_bn128_G2::wnaf_window_table.push_back(15);
    alt_bn128_G2::wnaf_window_table.push_back(39);
    alt_bn128_G2::wnaf_window_table.push_back(109);

    alt_bn128_G2::fixed_base_exp_window_table.resize(0);
    // window 1 is unbeaten in [-inf, 5.10]
    alt_bn128_G2::fixed_base_exp_window_table.push_back(1);
    // window 2 is unbeaten in [5.10, 10.43]
    alt_bn128_G2::fixed_base_exp_window_table.push_back(5);
    // window 3 is unbeaten in [10.43, 25.28]
    alt_bn128_G2::fixed_base_exp_window_table.push_back(10);
    // window 4 is unbeaten in [25.28, 59.00]
    alt_bn128_G2::fixed_base_exp_window_table.push_back(25);
    // window 5 is unbeaten in [59.00, 154.03]
    alt_bn128_G2::fixed_base_exp_window_table.push_back(59);
    // window 6 is unbeaten in [154.03, 334.25]
    alt_bn128_G2::fixed_base_exp_window_table.push_back(154);
    // window 7 is unbeaten in [334.25, 742.58]
    alt_bn128_G2::fixed_base_exp_window_table.push_back(334);
    // window 8 is unbeaten in [742.58, 2034.40]
    alt_bn128_G2::fixed_base_exp_window_table.push_back(743);
    // window 9 is unbeaten in [2034.40, 4987.56]
    alt_bn128_G2::fixed_base_exp_window_table.push_back(2034);
    // window 10 is unbeaten in [4987.56, 8888.27]
    alt_bn128_G2::fixed_base_exp_window_table.push_back(4988);
    // window 11 is unbeaten in [8888.27, 26271.13]
    alt_bn128_G2::fixed_base_exp_window_table.push_back(8888);
    // window 12 is unbeaten in [26271.13, 39768.20]
    alt_bn128_G2::fixed_base_exp_window_table.push_back(26271);
    // window 13 is unbeaten in [39768.20, 106275.75]
    alt_bn128_G2::fixed_base_exp_window_table.push_back(39768);
    // window 14 is unbeaten in [106275.75, 141703.40]
    alt_bn128_G2::fixed_base_exp_window_table.push_back(106276);
    // window 15 is unbeaten in [141703.40, 462422.97]
    alt_bn128_G2::fixed_base_exp_window_table.push_back(141703);
    // window 16 is unbeaten in [462422.97, 926871.84]
    alt_bn128_G2::fixed_base_exp_window_table.push_back(462423);
    // window 17 is unbeaten in [926871.84, 4873049.17]
    alt_bn128_G2::fixed_base_exp_window_table.push_back(926872);
    // window 18 is never the best
    alt_bn128_G2::fixed_base_exp_window_table.push_back(0);
    // window 19 is unbeaten in [4873049.17, 5706707.88]
    alt_bn128_G2::fixed_base_exp_window_table.push_back(4873049);
    // window 20 is unbeaten in [5706707.88, 31673814.95]
    alt_bn128_G2::fixed_base_exp_window_table.push_back(5706708);
    // window 21 is never the best
    alt_bn128_G2::fixed_base_exp_window_table.push_back(0);
    // window 22 is unbeaten in [31673814.95, inf]
    alt_bn128_G2::fixed_base_exp_window_table.push_back(31673815);

    /* pairing parameters */

    alt_bn128_ate_loop_count = bigint<alt_bn128_q_limbs>("29793968203157093288");
    alt_bn128_ate_is_loop_count_neg = false;
    alt_bn128_final_exponent = bigint<12*alt_bn128_q_limbs>("552484233613224096312617126783173147097382103762957654188882734314196910839907541213974502761540629817009608548654680343627701153829446747810907373256841551006201639677726139946029199968412598804882391702273019083653272047566316584365559776493027495458238373902875937659943504873220554161550525926302303331747463515644711876653177129578303191095900909191624817826566688241804408081892785725967931714097716709526092261278071952560171111444072049229123565057483750161460024353346284167282452756217662335528813519139808291170539072125381230815729071544861602750936964829313608137325426383735122175229541155376346436093930287402089517426973178917569713384748081827255472576937471496195752727188261435633271238710131736096299798168852925540549342330775279877006784354801422249722573783561685179618816480037695005515426162362431072245638324744480");
    alt_bn128_final_exponent_z = bigint<alt_bn128_q_limbs>("4965661367192848881");
    alt_bn128_final_exponent_is_z_neg = false;

}
} // libff
