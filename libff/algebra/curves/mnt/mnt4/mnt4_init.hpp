/** @file
 *****************************************************************************

 Declaration of interfaces for initializing MNT4.

 *****************************************************************************
 * @author     This file is part of libff, developed by SCIPR Lab
 *             and contributors (see AUTHORS).
 * @copyright  MIT license (see LICENSE file)
 *****************************************************************************/

#ifndef MNT4_INIT_HPP_
#define MNT4_INIT_HPP_

#include <libff/algebra/curves/public_params.hpp>
#include <libff/algebra/curves/mnt/mnt4/mnt4_fields.hpp>

namespace libff {

// parameters for twisted short Weierstrass curve E'/Fq2 : y^2 = x^3 + (a * twist^2) * x + (b * twist^3)
extern mnt4_Fq2 mnt4_twist;
extern mnt4_Fq2 mnt4_twist_coeff_a;
extern mnt4_Fq2 mnt4_twist_coeff_b;
extern mnt4_Fq mnt4_twist_mul_by_a_c0;
extern mnt4_Fq mnt4_twist_mul_by_a_c1;
extern mnt4_Fq mnt4_twist_mul_by_b_c0;
extern mnt4_Fq mnt4_twist_mul_by_b_c1;
extern mnt4_Fq mnt4_twist_mul_by_q_X;
extern mnt4_Fq mnt4_twist_mul_by_q_Y;

// parameters for pairing
extern bigint<mnt4_q_limbs> mnt4_ate_loop_count;
extern bool mnt4_ate_is_loop_count_neg;
extern bigint<4*mnt4_q_limbs> mnt4_final_exponent;
extern bigint<mnt4_q_limbs> mnt4_final_exponent_last_chunk_abs_of_w0;
extern bool mnt4_final_exponent_last_chunk_is_w0_neg;
extern bigint<mnt4_q_limbs> mnt4_final_exponent_last_chunk_w1;

void init_mnt4_params();

class mnt4_G1;
class mnt4_G2;

} // libff

#endif // MNT4_INIT_HPP_
