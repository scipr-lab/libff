/** @file
 *****************************************************************************
 * @author     This file is part of libff, developed by SCIPR Lab
 *             and contributors (see AUTHORS).
 * @copyright  MIT license (see LICENSE file)
 *****************************************************************************/

#ifndef BN128_FIELDS_HPP_
#define BN128_FIELDS_HPP_
#include "depends/ate-pairing/include/bn.h"
#include <libff/algebra/fields/prime_base/fp.hpp>

namespace libff {

const mp_size_t bn128_r_bitcount = 254;
const mp_size_t bn128_q_bitcount = 254;

const mp_size_t bn128_r_limbs = (bn128_r_bitcount+GMP_NUMB_BITS-1)/GMP_NUMB_BITS;
const mp_size_t bn128_q_limbs = (bn128_q_bitcount+GMP_NUMB_BITS-1)/GMP_NUMB_BITS;

extern bigint<bn128_r_limbs> bn128_modulus_r;
extern bigint<bn128_q_limbs> bn128_modulus_q;

typedef Fp_model<bn128_r_limbs, bn128_modulus_r> bn128_Fr;
typedef Fp_model<bn128_q_limbs, bn128_modulus_q> bn128_Fq;

void init_bn128_fields();

} // namespace libff
#endif // BN128_FIELDS_HPP_
