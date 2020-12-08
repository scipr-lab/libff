/** @file
 *****************************************************************************
 * @author     This file is part of libff, developed by SCIPR Lab
 *             and contributors (see AUTHORS).
 * @copyright  MIT license (see LICENSE file)
 *****************************************************************************/

#ifndef BN128_INIT_HPP_
#define BN128_INIT_HPP_
#include "ate-pairing/include/bn.h"

#include <libff/algebra/curves/public_params.hpp>
#include <libff/algebra/curves/bn128/bn128_fields.hpp>

namespace libff {

extern bn::Fp bn128_coeff_b;
extern std::size_t bn128_Fq_s;
extern bn::Fp bn128_Fq_nqr_to_t;
extern mie::Vuint bn128_Fq_t_minus_1_over_2;

extern bn::Fp2 bn128_twist_coeff_b;
extern std::size_t bn128_Fq2_s;
extern bn::Fp2 bn128_Fq2_nqr_to_t;
extern mie::Vuint bn128_Fq2_t_minus_1_over_2;

void init_bn128_params();

class bn128_G1;
class bn128_G2;
class bn128_GT;
typedef bn128_GT bn128_Fq12;

} // libff
#endif // BN128_INIT_HPP_
