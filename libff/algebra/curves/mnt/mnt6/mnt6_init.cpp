/** @file
 *****************************************************************************

 Implementation of interfaces for initializing MNT6.

 See mnt6_init.hpp .

 *****************************************************************************
 * @author     This file is part of libff, developed by SCIPR Lab
 *             and contributors (see AUTHORS).
 * @copyright  MIT license (see LICENSE file)
 *****************************************************************************/

#include <libff/algebra/curves/mnt/mnt6/mnt6_g1.hpp>
#include <libff/algebra/curves/mnt/mnt6/mnt6_g2.hpp>
#include <libff/algebra/curves/mnt/mnt6/mnt6_init.hpp>

namespace libff {

mnt6_Fq3 mnt6_twist;
mnt6_Fq3 mnt6_twist_coeff_a;
mnt6_Fq3 mnt6_twist_coeff_b;
mnt6_Fq mnt6_twist_mul_by_a_c0;
mnt6_Fq mnt6_twist_mul_by_a_c1;
mnt6_Fq mnt6_twist_mul_by_a_c2;
mnt6_Fq mnt6_twist_mul_by_b_c0;
mnt6_Fq mnt6_twist_mul_by_b_c1;
mnt6_Fq mnt6_twist_mul_by_b_c2;
mnt6_Fq mnt6_twist_mul_by_q_X;
mnt6_Fq mnt6_twist_mul_by_q_Y;

bigint<mnt6_q_limbs> mnt6_ate_loop_count;
bool mnt6_ate_is_loop_count_neg;
bigint<6*mnt6_q_limbs> mnt6_final_exponent;
bigint<mnt6_q_limbs> mnt6_final_exponent_last_chunk_abs_of_w0;
bool mnt6_final_exponent_last_chunk_is_w0_neg;
bigint<mnt6_q_limbs> mnt6_final_exponent_last_chunk_w1;

void init_mnt6_params()
{
    init_mnt6_fields();

    /* choice of short Weierstrass curve and its twist */
    mnt6_G1::coeff_a = mnt6_Fq("11");
    mnt6_G1::coeff_b = mnt6_Fq("106700080510851735677967319632585352256454251201367587890185989362936000262606668469523074");
    mnt6_twist = mnt6_Fq3(mnt6_Fq::zero(), mnt6_Fq::one(), mnt6_Fq::zero());
    mnt6_twist_coeff_a = mnt6_Fq3(mnt6_Fq::zero(), mnt6_Fq::zero(),
                                  mnt6_G1::coeff_a);
    mnt6_twist_coeff_b = mnt6_Fq3(mnt6_G1::coeff_b * mnt6_Fq3::non_residue,
                                  mnt6_Fq::zero(), mnt6_Fq::zero());
    mnt6_G2::twist = mnt6_twist;
    mnt6_G2::coeff_a = mnt6_twist_coeff_a;
    mnt6_G2::coeff_b = mnt6_twist_coeff_b;
    mnt6_twist_mul_by_a_c0 = mnt6_G1::coeff_a * mnt6_Fq3::non_residue;
    mnt6_twist_mul_by_a_c1 = mnt6_G1::coeff_a * mnt6_Fq3::non_residue;
    mnt6_twist_mul_by_a_c2 = mnt6_G1::coeff_a;
    mnt6_twist_mul_by_b_c0 = mnt6_G1::coeff_b * mnt6_Fq3::non_residue;
    mnt6_twist_mul_by_b_c1 = mnt6_G1::coeff_b * mnt6_Fq3::non_residue;
    mnt6_twist_mul_by_b_c2 = mnt6_G1::coeff_b * mnt6_Fq3::non_residue;
    mnt6_twist_mul_by_q_X = mnt6_Fq("4183387201740296620308398334599285547820769823264541783190415909159130177461911693276180");
    mnt6_twist_mul_by_q_Y = mnt6_Fq("475922286169261325753349249653048451545124878552823515553267735739164647307408490559963136");

    /* choice of group G1 */
    mnt6_G1::G1_zero = mnt6_G1(mnt6_Fq::zero(),
                               mnt6_Fq::one(),
                               mnt6_Fq::zero());
    mnt6_G1::G1_one = mnt6_G1(mnt6_Fq("336685752883082228109289846353937104185698209371404178342968838739115829740084426881123453"),
                              mnt6_Fq("402596290139780989709332707716568920777622032073762749862342374583908837063963736098549800"),
                              mnt6_Fq::one());
    mnt6_G1::initialized = true;

    // Cofactor
    mnt6_G1::h = bigint<mnt6_G1::h_limbs>("1");

    mnt6_G1::wnaf_window_table.resize(0);
    mnt6_G1::wnaf_window_table.push_back(11);
    mnt6_G1::wnaf_window_table.push_back(24);
    mnt6_G1::wnaf_window_table.push_back(60);
    mnt6_G1::wnaf_window_table.push_back(127);

    mnt6_G1::fixed_base_exp_window_table.resize(0);
    // window 1 is unbeaten in [-inf, 3.96]
    mnt6_G1::fixed_base_exp_window_table.push_back(1);
    // window 2 is unbeaten in [3.96, 9.67]
    mnt6_G1::fixed_base_exp_window_table.push_back(4);
    // window 3 is unbeaten in [9.67, 25.13]
    mnt6_G1::fixed_base_exp_window_table.push_back(10);
    // window 4 is unbeaten in [25.13, 60.31]
    mnt6_G1::fixed_base_exp_window_table.push_back(25);
    // window 5 is unbeaten in [60.31, 146.07]
    mnt6_G1::fixed_base_exp_window_table.push_back(60);
    // window 6 is unbeaten in [146.07, 350.09]
    mnt6_G1::fixed_base_exp_window_table.push_back(146);
    // window 7 is unbeaten in [350.09, 844.54]
    mnt6_G1::fixed_base_exp_window_table.push_back(350);
    // window 8 is unbeaten in [844.54, 1839.64]
    mnt6_G1::fixed_base_exp_window_table.push_back(845);
    // window 9 is unbeaten in [1839.64, 3904.26]
    mnt6_G1::fixed_base_exp_window_table.push_back(1840);
    // window 10 is unbeaten in [3904.26, 11309.42]
    mnt6_G1::fixed_base_exp_window_table.push_back(3904);
    // window 11 is unbeaten in [11309.42, 24015.57]
    mnt6_G1::fixed_base_exp_window_table.push_back(11309);
    // window 12 is unbeaten in [24015.57, 72288.57]
    mnt6_G1::fixed_base_exp_window_table.push_back(24016);
    // window 13 is unbeaten in [72288.57, 138413.22]
    mnt6_G1::fixed_base_exp_window_table.push_back(72289);
    // window 14 is unbeaten in [138413.22, 156390.30]
    mnt6_G1::fixed_base_exp_window_table.push_back(138413);
    // window 15 is unbeaten in [156390.30, 562560.50]
    mnt6_G1::fixed_base_exp_window_table.push_back(156390);
    // window 16 is unbeaten in [562560.50, 1036742.02]
    mnt6_G1::fixed_base_exp_window_table.push_back(562560);
    // window 17 is unbeaten in [1036742.02, 2053818.86]
    mnt6_G1::fixed_base_exp_window_table.push_back(1036742);
    // window 18 is unbeaten in [2053818.86, 4370223.95]
    mnt6_G1::fixed_base_exp_window_table.push_back(2053819);
    // window 19 is unbeaten in [4370223.95, 8215703.81]
    mnt6_G1::fixed_base_exp_window_table.push_back(4370224);
    // window 20 is unbeaten in [8215703.81, 42682375.43]
    mnt6_G1::fixed_base_exp_window_table.push_back(8215704);
    // window 21 is never the best
    mnt6_G1::fixed_base_exp_window_table.push_back(0);
    // window 22 is unbeaten in [42682375.43, inf]
    mnt6_G1::fixed_base_exp_window_table.push_back(42682375);

    /* choice of group G2 */
    mnt6_G2::G2_zero = mnt6_G2(mnt6_Fq3::zero(),
                               mnt6_Fq3::one(),
                               mnt6_Fq3::zero());
    mnt6_G2::G2_one = mnt6_G2(mnt6_Fq3(mnt6_Fq("421456435772811846256826561593908322288509115489119907560382401870203318738334702321297427"),
                                       mnt6_Fq("103072927438548502463527009961344915021167584706439945404959058962657261178393635706405114"),
                                       mnt6_Fq("143029172143731852627002926324735183809768363301149009204849580478324784395590388826052558")),
                              mnt6_Fq3(mnt6_Fq("464673596668689463130099227575639512541218133445388869383893594087634649237515554342751377"),
                                       mnt6_Fq("100642907501977375184575075967118071807821117960152743335603284583254620685343989304941678"),
                                       mnt6_Fq("123019855502969896026940545715841181300275180157288044663051565390506010149881373807142903")),
                              mnt6_Fq3::one());
    mnt6_G2::initialized = true;

    // Cofactor
    mnt6_G2::h = bigint<mnt6_G2::h_limbs>("226502022472576270196498690498308461791828762732602586162207535351960270082712694977333372361549082214519252261735048131889018501404377856786623430385820659037970876666767495659520");

    mnt6_G2::wnaf_window_table.resize(0);
    mnt6_G2::wnaf_window_table.push_back(5);
    mnt6_G2::wnaf_window_table.push_back(15);
    mnt6_G2::wnaf_window_table.push_back(39);
    mnt6_G2::wnaf_window_table.push_back(109);

    mnt6_G2::fixed_base_exp_window_table.resize(0);
    // window 1 is unbeaten in [-inf, 4.25]
    mnt6_G2::fixed_base_exp_window_table.push_back(1);
    // window 2 is unbeaten in [4.25, 10.22]
    mnt6_G2::fixed_base_exp_window_table.push_back(4);
    // window 3 is unbeaten in [10.22, 24.85]
    mnt6_G2::fixed_base_exp_window_table.push_back(10);
    // window 4 is unbeaten in [24.85, 60.06]
    mnt6_G2::fixed_base_exp_window_table.push_back(25);
    // window 5 is unbeaten in [60.06, 143.61]
    mnt6_G2::fixed_base_exp_window_table.push_back(60);
    // window 6 is unbeaten in [143.61, 345.66]
    mnt6_G2::fixed_base_exp_window_table.push_back(144);
    // window 7 is unbeaten in [345.66, 818.56]
    mnt6_G2::fixed_base_exp_window_table.push_back(346);
    // window 8 is unbeaten in [818.56, 1782.06]
    mnt6_G2::fixed_base_exp_window_table.push_back(819);
    // window 9 is unbeaten in [1782.06, 4002.45]
    mnt6_G2::fixed_base_exp_window_table.push_back(1782);
    // window 10 is unbeaten in [4002.45, 10870.18]
    mnt6_G2::fixed_base_exp_window_table.push_back(4002);
    // window 11 is unbeaten in [10870.18, 18022.51]
    mnt6_G2::fixed_base_exp_window_table.push_back(10870);
    // window 12 is unbeaten in [18022.51, 43160.74]
    mnt6_G2::fixed_base_exp_window_table.push_back(18023);
    // window 13 is unbeaten in [43160.74, 149743.32]
    mnt6_G2::fixed_base_exp_window_table.push_back(43161);
    // window 14 is never the best
    mnt6_G2::fixed_base_exp_window_table.push_back(0);
    // window 15 is unbeaten in [149743.32, 551844.13]
    mnt6_G2::fixed_base_exp_window_table.push_back(149743);
    // window 16 is unbeaten in [551844.13, 1041827.91]
    mnt6_G2::fixed_base_exp_window_table.push_back(551844);
    // window 17 is unbeaten in [1041827.91, 1977371.53]
    mnt6_G2::fixed_base_exp_window_table.push_back(1041828);
    // window 18 is unbeaten in [1977371.53, 3703619.51]
    mnt6_G2::fixed_base_exp_window_table.push_back(1977372);
    // window 19 is unbeaten in [3703619.51, 7057236.87]
    mnt6_G2::fixed_base_exp_window_table.push_back(3703620);
    // window 20 is unbeaten in [7057236.87, 38554491.67]
    mnt6_G2::fixed_base_exp_window_table.push_back(7057237);
    // window 21 is never the best
    mnt6_G2::fixed_base_exp_window_table.push_back(0);
    // window 22 is unbeaten in [38554491.67, inf]
    mnt6_G2::fixed_base_exp_window_table.push_back(38554492);

    /* pairing parameters */
    mnt6_ate_loop_count = bigint<mnt6_q_limbs>("689871209842287392837045615510547309923794944");
    mnt6_ate_is_loop_count_neg = true;
    mnt6_final_exponent = bigint<6*mnt6_q_limbs>("24416320138090509697890595414313438768353977489862543935904010715439066975957855922532159264213056712140358746422742237328406558352706591021642230618060502855451264045397444793186876199015256781648746888625527075466063075011307800862173764236311342105211681121426931616843635215852236649271569251468773714424208521977615548771268520882870120900360322044218806712027729351845307690474985502587527753847200130592058098363641559341826790559426614919168");
    mnt6_final_exponent_last_chunk_abs_of_w0 = bigint<mnt6_q_limbs>("689871209842287392837045615510547309923794944");
    mnt6_final_exponent_last_chunk_is_w0_neg = true;
    mnt6_final_exponent_last_chunk_w1 = bigint<mnt6_q_limbs>("1");
}

} // namespace libff
