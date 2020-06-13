/** @file
 *****************************************************************************

 Implementation of interfaces for initializing MNT4.

 See mnt4_init.hpp .

 *****************************************************************************
 * @author     This file is part of libff, developed by SCIPR Lab
 *             and contributors (see AUTHORS).
 * @copyright  MIT license (see LICENSE file)
 *****************************************************************************/

#include <libff/algebra/curves/mnt/mnt4/mnt4_g1.hpp>
#include <libff/algebra/curves/mnt/mnt4/mnt4_g2.hpp>
#include <libff/algebra/curves/mnt/mnt4/mnt4_init.hpp>

namespace libff {

mnt4_Fq2 mnt4_twist;
mnt4_Fq2 mnt4_twist_coeff_a;
mnt4_Fq2 mnt4_twist_coeff_b;
mnt4_Fq mnt4_twist_mul_by_a_c0;
mnt4_Fq mnt4_twist_mul_by_a_c1;
mnt4_Fq mnt4_twist_mul_by_b_c0;
mnt4_Fq mnt4_twist_mul_by_b_c1;
mnt4_Fq mnt4_twist_mul_by_q_X;
mnt4_Fq mnt4_twist_mul_by_q_Y;

bigint<mnt4_q_limbs> mnt4_ate_loop_count;
bool mnt4_ate_is_loop_count_neg;
bigint<4*mnt4_q_limbs> mnt4_final_exponent;
bigint<mnt4_q_limbs> mnt4_final_exponent_last_chunk_abs_of_w0;
bool mnt4_final_exponent_last_chunk_is_w0_neg;
bigint<mnt4_q_limbs> mnt4_final_exponent_last_chunk_w1;

void init_mnt4_params()
{
    init_mnt4_fields();

    /* choice of short Weierstrass curve and its twist */
    mnt4_G1::coeff_a = mnt4_Fq("2");
    mnt4_G1::coeff_b = mnt4_Fq("423894536526684178289416011533888240029318103673896002803341544124054745019340795360841685");
    mnt4_twist = mnt4_Fq2(mnt4_Fq::zero(), mnt4_Fq::one());
    mnt4_twist_coeff_a = mnt4_Fq2(mnt4_G1::coeff_a * mnt4_Fq2::non_residue, mnt4_Fq::zero());
    mnt4_twist_coeff_b = mnt4_Fq2(mnt4_Fq::zero(), mnt4_G1::coeff_b * mnt4_Fq2::non_residue);
    mnt4_G2::twist = mnt4_twist;
    mnt4_G2::coeff_a = mnt4_twist_coeff_a;
    mnt4_G2::coeff_b = mnt4_twist_coeff_b;
    mnt4_twist_mul_by_a_c0 = mnt4_G1::coeff_a * mnt4_Fq2::non_residue;
    mnt4_twist_mul_by_a_c1 = mnt4_G1::coeff_a * mnt4_Fq2::non_residue;
    mnt4_twist_mul_by_b_c0 = mnt4_G1::coeff_b * mnt4_Fq2::non_residue.squared();
    mnt4_twist_mul_by_b_c1 = mnt4_G1::coeff_b * mnt4_Fq2::non_residue;
    mnt4_twist_mul_by_q_X = mnt4_Fq("475922286169261325753349249653048451545124879242694725395555128576210262817955800483758080");
    mnt4_twist_mul_by_q_Y = mnt4_Fq("7684163245453501615621351552473337069301082060976805004625011694147890954040864167002308");

    /* choice of group G1 */
    mnt4_G1::G1_zero = mnt4_G1(mnt4_Fq::zero(),
                               mnt4_Fq::one(),
                               mnt4_Fq::zero());


    mnt4_G1::G1_one = mnt4_G1(mnt4_Fq("60760244141852568949126569781626075788424196370144486719385562369396875346601926534016838"),
                              mnt4_Fq("363732850702582978263902770815145784459747722357071843971107674179038674942891694705904306"),
                              mnt4_Fq::one());
    mnt4_G1::initialized = true;

    mnt4_G1::wnaf_window_table.resize(0);
    mnt4_G1::wnaf_window_table.push_back(11);
    mnt4_G1::wnaf_window_table.push_back(24);
    mnt4_G1::wnaf_window_table.push_back(60);
    mnt4_G1::wnaf_window_table.push_back(127);

    mnt4_G1::fixed_base_exp_window_table.resize(0);
    // window 1 is unbeaten in [-inf, 5.09]
    mnt4_G1::fixed_base_exp_window_table.push_back(1);
    // window 2 is unbeaten in [5.09, 9.64]
    mnt4_G1::fixed_base_exp_window_table.push_back(5);
    // window 3 is unbeaten in [9.64, 24.79]
    mnt4_G1::fixed_base_exp_window_table.push_back(10);
    // window 4 is unbeaten in [24.79, 60.29]
    mnt4_G1::fixed_base_exp_window_table.push_back(25);
    // window 5 is unbeaten in [60.29, 144.37]
    mnt4_G1::fixed_base_exp_window_table.push_back(60);
    // window 6 is unbeaten in [144.37, 344.90]
    mnt4_G1::fixed_base_exp_window_table.push_back(144);
    // window 7 is unbeaten in [344.90, 855.00]
    mnt4_G1::fixed_base_exp_window_table.push_back(345);
    // window 8 is unbeaten in [855.00, 1804.62]
    mnt4_G1::fixed_base_exp_window_table.push_back(855);
    // window 9 is unbeaten in [1804.62, 3912.30]
    mnt4_G1::fixed_base_exp_window_table.push_back(1805);
    // window 10 is unbeaten in [3912.30, 11264.50]
    mnt4_G1::fixed_base_exp_window_table.push_back(3912);
    // window 11 is unbeaten in [11264.50, 27897.51]
    mnt4_G1::fixed_base_exp_window_table.push_back(11265);
    // window 12 is unbeaten in [27897.51, 57596.79]
    mnt4_G1::fixed_base_exp_window_table.push_back(27898);
    // window 13 is unbeaten in [57596.79, 145298.71]
    mnt4_G1::fixed_base_exp_window_table.push_back(57597);
    // window 14 is unbeaten in [145298.71, 157204.59]
    mnt4_G1::fixed_base_exp_window_table.push_back(145299);
    // window 15 is unbeaten in [157204.59, 601600.62]
    mnt4_G1::fixed_base_exp_window_table.push_back(157205);
    // window 16 is unbeaten in [601600.62, 1107377.25]
    mnt4_G1::fixed_base_exp_window_table.push_back(601601);
    // window 17 is unbeaten in [1107377.25, 1789646.95]
    mnt4_G1::fixed_base_exp_window_table.push_back(1107377);
    // window 18 is unbeaten in [1789646.95, 4392626.92]
    mnt4_G1::fixed_base_exp_window_table.push_back(1789647);
    // window 19 is unbeaten in [4392626.92, 8221210.60]
    mnt4_G1::fixed_base_exp_window_table.push_back(4392627);
    // window 20 is unbeaten in [8221210.60, 42363731.19]
    mnt4_G1::fixed_base_exp_window_table.push_back(8221211);
    // window 21 is never the best
    mnt4_G1::fixed_base_exp_window_table.push_back(0);
    // window 22 is unbeaten in [42363731.19, inf]
    mnt4_G1::fixed_base_exp_window_table.push_back(42363731);

    /* choice of group G2 */
    mnt4_G2::G2_zero = mnt4_G2(mnt4_Fq2::zero(),
                               mnt4_Fq2::one(),
                               mnt4_Fq2::zero());

    mnt4_G2::G2_one = mnt4_G2(mnt4_Fq2(mnt4_Fq("438374926219350099854919100077809681842783509163790991847867546339851681564223481322252708"),
                                       mnt4_Fq("37620953615500480110935514360923278605464476459712393277679280819942849043649216370485641")),
                              mnt4_Fq2(mnt4_Fq("37437409008528968268352521034936931842973546441370663118543015118291998305624025037512482"),
                                       mnt4_Fq("424621479598893882672393190337420680597584695892317197646113820787463109735345923009077489")),
                              mnt4_Fq2::one());
    mnt4_G2::initialized = true;

    mnt4_G2::wnaf_window_table.resize(0);
    mnt4_G2::wnaf_window_table.push_back(5);
    mnt4_G2::wnaf_window_table.push_back(15);
    mnt4_G2::wnaf_window_table.push_back(39);
    mnt4_G2::wnaf_window_table.push_back(109);

    mnt4_G2::fixed_base_exp_window_table.resize(0);
    // window 1 is unbeaten in [-inf, 4.17]
    mnt4_G2::fixed_base_exp_window_table.push_back(1);
    // window 2 is unbeaten in [4.17, 10.12]
    mnt4_G2::fixed_base_exp_window_table.push_back(4);
    // window 3 is unbeaten in [10.12, 24.65]
    mnt4_G2::fixed_base_exp_window_table.push_back(10);
    // window 4 is unbeaten in [24.65, 60.03]
    mnt4_G2::fixed_base_exp_window_table.push_back(25);
    // window 5 is unbeaten in [60.03, 143.16]
    mnt4_G2::fixed_base_exp_window_table.push_back(60);
    // window 6 is unbeaten in [143.16, 344.73]
    mnt4_G2::fixed_base_exp_window_table.push_back(143);
    // window 7 is unbeaten in [344.73, 821.24]
    mnt4_G2::fixed_base_exp_window_table.push_back(345);
    // window 8 is unbeaten in [821.24, 1793.92]
    mnt4_G2::fixed_base_exp_window_table.push_back(821);
    // window 9 is unbeaten in [1793.92, 3919.59]
    mnt4_G2::fixed_base_exp_window_table.push_back(1794);
    // window 10 is unbeaten in [3919.59, 11301.46]
    mnt4_G2::fixed_base_exp_window_table.push_back(3920);
    // window 11 is unbeaten in [11301.46, 18960.09]
    mnt4_G2::fixed_base_exp_window_table.push_back(11301);
    // window 12 is unbeaten in [18960.09, 44198.62]
    mnt4_G2::fixed_base_exp_window_table.push_back(18960);
    // window 13 is unbeaten in [44198.62, 150799.57]
    mnt4_G2::fixed_base_exp_window_table.push_back(44199);
    // window 14 is never the best
    mnt4_G2::fixed_base_exp_window_table.push_back(0);
    // window 15 is unbeaten in [150799.57, 548694.81]
    mnt4_G2::fixed_base_exp_window_table.push_back(150800);
    // window 16 is unbeaten in [548694.81, 1051769.08]
    mnt4_G2::fixed_base_exp_window_table.push_back(548695);
    // window 17 is unbeaten in [1051769.08, 2023925.59]
    mnt4_G2::fixed_base_exp_window_table.push_back(1051769);
    // window 18 is unbeaten in [2023925.59, 3787108.68]
    mnt4_G2::fixed_base_exp_window_table.push_back(2023926);
    // window 19 is unbeaten in [3787108.68, 7107480.30]
    mnt4_G2::fixed_base_exp_window_table.push_back(3787109);
    // window 20 is unbeaten in [7107480.30, 38760027.14]
    mnt4_G2::fixed_base_exp_window_table.push_back(7107480);
    // window 21 is never the best
    mnt4_G2::fixed_base_exp_window_table.push_back(0);
    // window 22 is unbeaten in [38760027.14, inf]
    mnt4_G2::fixed_base_exp_window_table.push_back(38760027);

    /* pairing parameters */
    mnt4_ate_loop_count = bigint<mnt4_q_limbs>("689871209842287392837045615510547309923794944");
    mnt4_ate_is_loop_count_neg = false;
    mnt4_final_exponent = bigint<4*mnt4_q_limbs>("107797360357109903430794490309592072278927783803031854357910908121903439838772861497177116410825586743089760869945394610511917274977971559062689561855016270594656570874331111995170645233717143416875749097203441437192367065467706065411650403684877366879441766585988546560");
    mnt4_final_exponent_last_chunk_abs_of_w0 = bigint<mnt4_q_limbs>("689871209842287392837045615510547309923794945");
    mnt4_final_exponent_last_chunk_is_w0_neg = false;
    mnt4_final_exponent_last_chunk_w1 = bigint<mnt4_q_limbs>("1");
}

} // libff
