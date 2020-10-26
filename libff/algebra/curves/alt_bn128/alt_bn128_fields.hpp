/** @file
 *****************************************************************************
 * @author     This file is part of libff, developed by SCIPR Lab
 *             and contributors (see AUTHORS).
 * @copyright  MIT license (see LICENSE file)
 *****************************************************************************/

#ifndef ALT_BN128_FIELDS_HPP_
#define ALT_BN128_FIELDS_HPP_
#include <libff/algebra/fields/prime_base/fp.hpp>
#include <libff/algebra/fields/prime_extension/fp12_2over3over2.hpp>
#include <libff/algebra/fields/prime_extension/fp2.hpp>
#include <libff/algebra/fields/prime_extension/fp6_3over2.hpp>

namespace libff {

const mp_size_t alt_bn128_r_bitcount = 254;
const mp_size_t alt_bn128_q_bitcount = 254;

const mp_size_t alt_bn128_r_limbs = (alt_bn128_r_bitcount+GMP_NUMB_BITS-1)/GMP_NUMB_BITS;
const mp_size_t alt_bn128_q_limbs = (alt_bn128_q_bitcount+GMP_NUMB_BITS-1)/GMP_NUMB_BITS;

extern bigint<alt_bn128_r_limbs> alt_bn128_modulus_r;
extern bigint<alt_bn128_q_limbs> alt_bn128_modulus_q;

typedef Fp_model<alt_bn128_r_limbs, alt_bn128_modulus_r> alt_bn128_Fr;
typedef Fp_model<alt_bn128_q_limbs, alt_bn128_modulus_q> alt_bn128_Fq;
typedef Fp2_model<alt_bn128_q_limbs, alt_bn128_modulus_q> alt_bn128_Fq2;
typedef Fp6_3over2_model<alt_bn128_q_limbs, alt_bn128_modulus_q> alt_bn128_Fq6;
typedef Fp12_2over3over2_model<alt_bn128_q_limbs, alt_bn128_modulus_q> alt_bn128_Fq12;
typedef alt_bn128_Fq12 alt_bn128_GT;

void init_alt_bn128_fields();

} // libff
#endif // ALT_BN128_FIELDS_HPP_
