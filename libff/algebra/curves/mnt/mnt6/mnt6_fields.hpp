/** @file
 *****************************************************************************

 Declaration of interfaces for initializing MNT6.

 *****************************************************************************
 * @author     This file is part of libff, developed by SCIPR Lab
 *             and contributors (see AUTHORS).
 * @copyright  MIT license (see LICENSE file)
 *****************************************************************************/

#ifndef MNT6_FIELDS_HPP_
#define MNT6_FIELDS_HPP_

#include <libff/algebra/curves/mnt/mnt46_common.hpp>
#include <libff/algebra/curves/public_params.hpp>
#include <libff/algebra/fields/prime_base/fp.hpp>
#include <libff/algebra/fields/prime_extension/fp3.hpp>
#include <libff/algebra/fields/prime_extension/fp6_2over3.hpp>

namespace libff {

#define mnt6_modulus_r mnt46_modulus_B
#define mnt6_modulus_q mnt46_modulus_A

const mp_size_t mnt6_r_bitcount = mnt46_B_bitcount;
const mp_size_t mnt6_q_bitcount = mnt46_A_bitcount;

const mp_size_t mnt6_r_limbs = mnt46_B_limbs;
const mp_size_t mnt6_q_limbs = mnt46_A_limbs;

extern bigint<mnt6_r_limbs> mnt6_modulus_r;
extern bigint<mnt6_q_limbs> mnt6_modulus_q;

typedef Fp_model<mnt6_r_limbs, mnt6_modulus_r> mnt6_Fr;
typedef Fp_model<mnt6_q_limbs, mnt6_modulus_q> mnt6_Fq;
typedef Fp3_model<mnt6_q_limbs, mnt6_modulus_q> mnt6_Fq3;
typedef Fp6_2over3_model<mnt6_q_limbs, mnt6_modulus_q> mnt6_Fq6;
typedef mnt6_Fq6 mnt6_GT;

void init_mnt6_fields();

} // namespace libff

#endif // MNT6_FIELDS_HPP_
