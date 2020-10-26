/** @file
 *****************************************************************************
 * @author     This file is part of libff, developed by SCIPR Lab
 *             and contributors (see AUTHORS).
 * @copyright  MIT license (see LICENSE file)
 *****************************************************************************/

#ifndef EDWARDS_INIT_HPP_
#define EDWARDS_INIT_HPP_
#include <libff/algebra/curves/public_params.hpp>
#include <libff/algebra/curves/edwards/edwards_fields.hpp>

namespace libff {

// parameters for Edwards curve E_{1,d}(F_q)
extern edwards_Fq edwards_coeff_a;
extern edwards_Fq edwards_coeff_d;
// parameters for twisted Edwards curve E_{a',d'}(F_q^3)
extern edwards_Fq3 edwards_twist;
extern edwards_Fq3 edwards_twist_coeff_a;
extern edwards_Fq3 edwards_twist_coeff_d;
extern edwards_Fq edwards_twist_mul_by_a_c0;
extern edwards_Fq edwards_twist_mul_by_a_c1;
extern edwards_Fq edwards_twist_mul_by_a_c2;
extern edwards_Fq edwards_twist_mul_by_d_c0;
extern edwards_Fq edwards_twist_mul_by_d_c1;
extern edwards_Fq edwards_twist_mul_by_d_c2;
extern edwards_Fq edwards_twist_mul_by_q_Y;
extern edwards_Fq edwards_twist_mul_by_q_Z;

// parameters for pairing
extern bigint<edwards_q_limbs> edwards_ate_loop_count;
extern bigint<6*edwards_q_limbs> edwards_final_exponent;
extern bigint<edwards_q_limbs> edwards_final_exponent_last_chunk_abs_of_w0;
extern bool edwards_final_exponent_last_chunk_is_w0_neg;
extern bigint<edwards_q_limbs> edwards_final_exponent_last_chunk_w1;

void init_edwards_params();

class edwards_G1;
class edwards_G2;

} // libff
#endif // EDWARDS_INIT_HPP_
