#ifndef BLS12_381_INIT_HPP_
#define BLS12_381_INIT_HPP_
#include <libff/algebra/curves/public_params.hpp>
#include <libff/algebra/curves/bls12_381/bls12_381_fields.hpp>

namespace libff {

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
