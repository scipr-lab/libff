/** @file
 *****************************************************************************

 Declaration of interfaces for initializing MNT4.

 *****************************************************************************
 * @author     This file is part of libff, developed by SCIPR Lab
 *             and contributors (see AUTHORS).
 * @copyright  MIT license (see LICENSE file)
 *****************************************************************************/

#ifndef MNT4_FIELDS_HPP_
#define MNT4_FIELDS_HPP_

#include <libff/algebra/curves/mnt/mnt46_common.hpp>
#include <libff/algebra/fields/prime_base/fp.hpp>
#include <libff/algebra/fields/prime_extension/fp2.hpp>
#include <libff/algebra/fields/prime_extension/fp4.hpp>

namespace libff {

#define mnt4_modulus_r mnt46_modulus_A
#define mnt4_modulus_q mnt46_modulus_B

const mp_size_t mnt4_r_bitcount = mnt46_A_bitcount;
const mp_size_t mnt4_q_bitcount = mnt46_B_bitcount;

const mp_size_t mnt4_r_limbs = mnt46_A_limbs;
const mp_size_t mnt4_q_limbs = mnt46_B_limbs;

extern bigint<mnt4_r_limbs> mnt4_modulus_r;
extern bigint<mnt4_q_limbs> mnt4_modulus_q;

typedef Fp_model<mnt4_r_limbs, mnt4_modulus_r> mnt4_Fr;
typedef Fp_model<mnt4_q_limbs, mnt4_modulus_q> mnt4_Fq;
typedef Fp2_model<mnt4_q_limbs, mnt4_modulus_q> mnt4_Fq2;
typedef Fp4_model<mnt4_q_limbs, mnt4_modulus_q> mnt4_Fq4;
typedef mnt4_Fq4 mnt4_GT;

void init_mnt4_fields();

} // libff

#endif // MNT4_FIELDS_HPP_
