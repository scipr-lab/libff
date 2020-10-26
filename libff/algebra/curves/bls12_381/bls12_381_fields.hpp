/** @file
 *****************************************************************************
 * @author     This file is part of libff, developed by SCIPR Lab
 *             and contributors (see AUTHORS).
 * @copyright  MIT license (see LICENSE file)
 *****************************************************************************/

#ifndef BLS12_381_FIELDS_HPP_
#define BLS12_381_FIELDS_HPP_
#include <libff/algebra/fields/prime_base/fp.hpp>
#include <libff/algebra/fields/prime_extension/fp2.hpp>
#include <libff/algebra/fields/prime_extension/fp6_3over2.hpp>
#include <libff/algebra/fields/prime_extension/fp12_2over3over2.hpp>

namespace libff {

const mp_size_t bls12_381_r_bitcount = 255;
const mp_size_t bls12_381_q_bitcount = 381;

const mp_size_t bls12_381_r_limbs = (bls12_381_r_bitcount+GMP_NUMB_BITS-1)/GMP_NUMB_BITS;
const mp_size_t bls12_381_q_limbs = (bls12_381_q_bitcount+GMP_NUMB_BITS-1)/GMP_NUMB_BITS;

extern bigint<bls12_381_r_limbs> bls12_381_modulus_r;
extern bigint<bls12_381_q_limbs> bls12_381_modulus_q;

typedef Fp_model<bls12_381_r_limbs, bls12_381_modulus_r> bls12_381_Fr;
typedef Fp_model<bls12_381_q_limbs, bls12_381_modulus_q> bls12_381_Fq;
typedef Fp2_model<bls12_381_q_limbs, bls12_381_modulus_q> bls12_381_Fq2;
typedef Fp6_3over2_model<bls12_381_q_limbs, bls12_381_modulus_q> bls12_381_Fq6;
typedef Fp12_2over3over2_model<bls12_381_q_limbs, bls12_381_modulus_q> bls12_381_Fq12;
typedef bls12_381_Fq12 bls12_381_GT;

void init_bls12_381_fields();

} // libff
#endif // BLS12_381_FIELDS_HPP_
