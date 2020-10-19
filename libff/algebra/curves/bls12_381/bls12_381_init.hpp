#ifndef BLS12_381_INIT_HPP_
#define BLS12_381_INIT_HPP_
#include <libff/algebra/curves/public_params.hpp>
#include <libff/algebra/fields/fp.hpp>
#include <libff/algebra/fields/fp12_2over3over2.hpp>
#include <libff/algebra/fields/fp2.hpp>
#include <libff/algebra/fields/fp6_3over2.hpp>

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

// parameters for the curve E/Fq : y^2 = x^3 + b
extern bls12_381_Fq bls12_381_coeff_b;
// parameters for the twisted curve E'/Fq2 : y^2 = x^3 + b/xi
extern bls12_381_Fq2 bls12_381_twist;
extern bls12_381_Fq2 bls12_381_twist_coeff_b;
extern bls12_381_Fq bls12_381_twist_mul_by_b_c0;
extern bls12_381_Fq bls12_381_twist_mul_by_b_c1;
extern bls12_381_Fq2 bls12_381_twist_mul_by_q_X;
extern bls12_381_Fq2 bls12_381_twist_mul_by_q_Y;

// parameters for pairing
extern bigint<bls12_381_q_limbs> bls12_381_ate_loop_count;
extern bool bls12_381_ate_is_loop_count_neg;
extern bigint<12*bls12_381_q_limbs> bls12_381_final_exponent;
extern bigint<bls12_381_q_limbs> bls12_381_final_exponent_z;
extern bool bls12_381_final_exponent_is_z_neg;

void init_bls12_381_params();

class bls12_381_G1;
class bls12_381_G2;

} // libff
#endif // BLS12_381_INIT_HPP_
