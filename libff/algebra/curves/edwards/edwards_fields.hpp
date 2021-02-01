/** @file
 *****************************************************************************
 * @author     This file is part of libff, developed by SCIPR Lab
 *             and contributors (see AUTHORS).
 * @copyright  MIT license (see LICENSE file)
 *****************************************************************************/

#ifndef EDWARDS_FIELDS_HPP_
#define EDWARDS_FIELDS_HPP_
#include <libff/algebra/fields/prime_base/fp.hpp>
#include <libff/algebra/fields/prime_extension/fp3.hpp>
#include <libff/algebra/fields/prime_extension/fp6_2over3.hpp>

namespace libff {

const mp_size_t edwards_r_bitcount = 181;
const mp_size_t edwards_q_bitcount = 183;

const mp_size_t edwards_r_limbs = (edwards_r_bitcount+GMP_NUMB_BITS-1)/GMP_NUMB_BITS;
const mp_size_t edwards_q_limbs = (edwards_q_bitcount+GMP_NUMB_BITS-1)/GMP_NUMB_BITS;

extern bigint<edwards_r_limbs> edwards_modulus_r;
extern bigint<edwards_q_limbs> edwards_modulus_q;

typedef Fp_model<edwards_r_limbs, edwards_modulus_r> edwards_Fr;
typedef Fp_model<edwards_q_limbs, edwards_modulus_q> edwards_Fq;
typedef Fp3_model<edwards_q_limbs, edwards_modulus_q> edwards_Fq3;
typedef Fp6_2over3_model<edwards_q_limbs, edwards_modulus_q> edwards_Fq6;
typedef edwards_Fq6 edwards_GT;

void init_edwards_fields();

} // namespace libff
#endif // EDWARDS_FIELDS_HPP_
