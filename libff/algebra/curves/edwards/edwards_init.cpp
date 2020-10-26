/** @file
 *****************************************************************************
 * @author     This file is part of libff, developed by SCIPR Lab
 *             and contributors (see AUTHORS).
 * @copyright  MIT license (see LICENSE file)
 *****************************************************************************/

#include <libff/algebra/curves/edwards/edwards_g1.hpp>
#include <libff/algebra/curves/edwards/edwards_g2.hpp>
#include <libff/algebra/curves/edwards/edwards_init.hpp>

namespace libff {

edwards_Fq edwards_coeff_a;
edwards_Fq edwards_coeff_d;
edwards_Fq3 edwards_twist;
edwards_Fq3 edwards_twist_coeff_a;
edwards_Fq3 edwards_twist_coeff_d;
edwards_Fq edwards_twist_mul_by_a_c0;
edwards_Fq edwards_twist_mul_by_a_c1;
edwards_Fq edwards_twist_mul_by_a_c2;
edwards_Fq edwards_twist_mul_by_d_c0;
edwards_Fq edwards_twist_mul_by_d_c1;
edwards_Fq edwards_twist_mul_by_d_c2;
edwards_Fq edwards_twist_mul_by_q_Y;
edwards_Fq edwards_twist_mul_by_q_Z;

bigint<edwards_q_limbs> edwards_ate_loop_count;
bigint<6*edwards_q_limbs> edwards_final_exponent;
bigint<edwards_q_limbs> edwards_final_exponent_last_chunk_abs_of_w0;
bool edwards_final_exponent_last_chunk_is_w0_neg;
bigint<edwards_q_limbs> edwards_final_exponent_last_chunk_w1;

void init_edwards_params()
{
    init_edwards_fields();

    /* choice of Edwards curve and its twist */

    edwards_coeff_a = edwards_Fq::one();
    edwards_coeff_d = edwards_Fq("600581931845324488256649384912508268813600056237543024");
    edwards_twist = edwards_Fq3(edwards_Fq::zero(), edwards_Fq::one(), edwards_Fq::zero());
    edwards_twist_coeff_a = edwards_coeff_a * edwards_twist;
    edwards_twist_coeff_d = edwards_coeff_d * edwards_twist;
    edwards_twist_mul_by_a_c0 = edwards_coeff_a * edwards_Fq3::non_residue;
    edwards_twist_mul_by_a_c1 = edwards_coeff_a;
    edwards_twist_mul_by_a_c2 = edwards_coeff_a;
    edwards_twist_mul_by_d_c0 = edwards_coeff_d * edwards_Fq3::non_residue;
    edwards_twist_mul_by_d_c1 = edwards_coeff_d;
    edwards_twist_mul_by_d_c2 = edwards_coeff_d;
    edwards_twist_mul_by_q_Y = edwards_Fq("1073752683758513276629212192812154536507607213288832062");
    edwards_twist_mul_by_q_Z = edwards_Fq("1073752683758513276629212192812154536507607213288832062");

    /* choice of group G1 */
    edwards_G1::G1_zero = edwards_G1(edwards_Fq::zero(),
                                     edwards_Fq::one());
    edwards_G1::G1_one = edwards_G1(edwards_Fq("3713709671941291996998665608188072510389821008693530490"),
                                    edwards_Fq("4869953702976555123067178261685365085639705297852816679"));
    edwards_G1::initialized = true;

    edwards_G1::wnaf_window_table.resize(0);
    edwards_G1::wnaf_window_table.push_back(9);
    edwards_G1::wnaf_window_table.push_back(14);
    edwards_G1::wnaf_window_table.push_back(24);
    edwards_G1::wnaf_window_table.push_back(117);

    edwards_G1::fixed_base_exp_window_table.resize(0);
    // window 1 is unbeaten in [-inf, 4.10]
    edwards_G1::fixed_base_exp_window_table.push_back(1);
    // window 2 is unbeaten in [4.10, 9.69]
    edwards_G1::fixed_base_exp_window_table.push_back(4);
    // window 3 is unbeaten in [9.69, 25.21]
    edwards_G1::fixed_base_exp_window_table.push_back(10);
    // window 4 is unbeaten in [25.21, 60.00]
    edwards_G1::fixed_base_exp_window_table.push_back(25);
    // window 5 is unbeaten in [60.00, 149.33]
    edwards_G1::fixed_base_exp_window_table.push_back(60);
    // window 6 is unbeaten in [149.33, 369.61]
    edwards_G1::fixed_base_exp_window_table.push_back(149);
    // window 7 is unbeaten in [369.61, 849.07]
    edwards_G1::fixed_base_exp_window_table.push_back(370);
    // window 8 is unbeaten in [849.07, 1764.94]
    edwards_G1::fixed_base_exp_window_table.push_back(849);
    // window 9 is unbeaten in [1764.94, 4429.59]
    edwards_G1::fixed_base_exp_window_table.push_back(1765);
    // window 10 is unbeaten in [4429.59, 13388.78]
    edwards_G1::fixed_base_exp_window_table.push_back(4430);
    // window 11 is unbeaten in [13388.78, 15368.00]
    edwards_G1::fixed_base_exp_window_table.push_back(13389);
    // window 12 is unbeaten in [15368.00, 74912.07]
    edwards_G1::fixed_base_exp_window_table.push_back(15368);
    // window 13 is unbeaten in [74912.07, 438107.20]
    edwards_G1::fixed_base_exp_window_table.push_back(74912);
    // window 14 is never the best
    edwards_G1::fixed_base_exp_window_table.push_back(0);
    // window 15 is unbeaten in [438107.20, 1045626.18]
    edwards_G1::fixed_base_exp_window_table.push_back(438107);
    // window 16 is never the best
    edwards_G1::fixed_base_exp_window_table.push_back(0);
    // window 17 is unbeaten in [1045626.18, 1577434.48]
    edwards_G1::fixed_base_exp_window_table.push_back(1045626);
    // window 18 is unbeaten in [1577434.48, 17350594.23]
    edwards_G1::fixed_base_exp_window_table.push_back(1577434);
    // window 19 is never the best
    edwards_G1::fixed_base_exp_window_table.push_back(0);
    // window 20 is never the best
    edwards_G1::fixed_base_exp_window_table.push_back(0);
    // window 21 is unbeaten in [17350594.23, inf]
    edwards_G1::fixed_base_exp_window_table.push_back(17350594);
    // window 22 is never the best
    edwards_G1::fixed_base_exp_window_table.push_back(0);

    /* choice of group G2 */

    edwards_G2::G2_zero = edwards_G2(edwards_Fq3::zero(),
                                     edwards_Fq3::one());
    edwards_G2::G2_one = edwards_G2(edwards_Fq3(edwards_Fq("4531683359223370252210990718516622098304721701253228128"),
                                                edwards_Fq("5339624155305731263217400504407647531329993548123477368"),
                                                edwards_Fq("3964037981777308726208525982198654699800283729988686552")),
                                    edwards_Fq3(edwards_Fq("364634864866983740775341816274081071386963546650700569"),
                                                edwards_Fq("3264380230116139014996291397901297105159834497864380415"),
                                                edwards_Fq("3504781284999684163274269077749440837914479176282903747")));
    edwards_G2::initialized = true;

    edwards_G2::wnaf_window_table.resize(0);
    edwards_G2::wnaf_window_table.push_back(6);
    edwards_G2::wnaf_window_table.push_back(12);
    edwards_G2::wnaf_window_table.push_back(42);
    edwards_G2::wnaf_window_table.push_back(97);

    edwards_G2::fixed_base_exp_window_table.resize(0);
    // window 1 is unbeaten in [-inf, 4.74]
    edwards_G2::fixed_base_exp_window_table.push_back(1);
    // window 2 is unbeaten in [4.74, 10.67]
    edwards_G2::fixed_base_exp_window_table.push_back(5);
    // window 3 is unbeaten in [10.67, 25.53]
    edwards_G2::fixed_base_exp_window_table.push_back(11);
    // window 4 is unbeaten in [25.53, 60.67]
    edwards_G2::fixed_base_exp_window_table.push_back(26);
    // window 5 is unbeaten in [60.67, 145.77]
    edwards_G2::fixed_base_exp_window_table.push_back(61);
    // window 6 is unbeaten in [145.77, 356.76]
    edwards_G2::fixed_base_exp_window_table.push_back(146);
    // window 7 is unbeaten in [356.76, 823.08]
    edwards_G2::fixed_base_exp_window_table.push_back(357);
    // window 8 is unbeaten in [823.08, 1589.45]
    edwards_G2::fixed_base_exp_window_table.push_back(823);
    // window 9 is unbeaten in [1589.45, 4135.70]
    edwards_G2::fixed_base_exp_window_table.push_back(1589);
    // window 10 is unbeaten in [4135.70, 14297.74]
    edwards_G2::fixed_base_exp_window_table.push_back(4136);
    // window 11 is unbeaten in [14297.74, 16744.85]
    edwards_G2::fixed_base_exp_window_table.push_back(14298);
    // window 12 is unbeaten in [16744.85, 51768.98]
    edwards_G2::fixed_base_exp_window_table.push_back(16745);
    // window 13 is unbeaten in [51768.98, 99811.01]
    edwards_G2::fixed_base_exp_window_table.push_back(51769);
    // window 14 is unbeaten in [99811.01, 193306.72]
    edwards_G2::fixed_base_exp_window_table.push_back(99811);
    // window 15 is unbeaten in [193306.72, 907184.68]
    edwards_G2::fixed_base_exp_window_table.push_back(193307);
    // window 16 is never the best
    edwards_G2::fixed_base_exp_window_table.push_back(0);
    // window 17 is unbeaten in [907184.68, 1389682.59]
    edwards_G2::fixed_base_exp_window_table.push_back(907185);
    // window 18 is unbeaten in [1389682.59, 6752695.74]
    edwards_G2::fixed_base_exp_window_table.push_back(1389683);
    // window 19 is never the best
    edwards_G2::fixed_base_exp_window_table.push_back(0);
    // window 20 is unbeaten in [6752695.74, 193642894.51]
    edwards_G2::fixed_base_exp_window_table.push_back(6752696);
    // window 21 is unbeaten in [193642894.51, 226760202.29]
    edwards_G2::fixed_base_exp_window_table.push_back(193642895);
    // window 22 is unbeaten in [226760202.29, inf]
    edwards_G2::fixed_base_exp_window_table.push_back(226760202);

    /* pairing parameters */

    edwards_ate_loop_count = bigint<edwards_q_limbs>("4492509698523932320491110403");
    edwards_final_exponent = bigint<6*edwards_q_limbs>("36943107177961694649618797346446870138748651578611748415128207429491593976636391130175425245705674550269561361208979548749447898941828686017765730419416875539615941651269793928962468899856083169227457503942470721108165443528513330156264699608120624990672333642644221591552000");
    edwards_final_exponent_last_chunk_abs_of_w0 = bigint<edwards_q_limbs>("17970038794095729281964441603");
    edwards_final_exponent_last_chunk_is_w0_neg = true;
    edwards_final_exponent_last_chunk_w1 = bigint<edwards_q_limbs>("4");

}
} // libff
