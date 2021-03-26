/** @file
 *****************************************************************************
 * @author     This file is part of libff, developed by SCIPR Lab
 *             and contributors (see AUTHORS).
 * @copyright  MIT license (see LICENSE file)
 *****************************************************************************/

#ifndef BN128_G2_HPP_
#define BN128_G2_HPP_
#include <iostream>
#include <vector>

#include "depends/ate-pairing/include/bn.h"

#include <libff/algebra/curves/bn128/bn128_init.hpp>
#include <libff/algebra/curves/curve_utils.hpp>

namespace libff {

class bn128_G2;
std::ostream& operator<<(std::ostream &, const bn128_G2&);
std::istream& operator>>(std::istream &, bn128_G2&);

class bn128_G2 {
private:
    static bn::Fp2 sqrt(const bn::Fp2 &el);
public:
#ifdef PROFILE_OP_COUNTS
    static long long add_cnt;
    static long long dbl_cnt;
#endif
    static std::vector<std::size_t> wnaf_window_table;
    static std::vector<std::size_t> fixed_base_exp_window_table;
    static bn128_G2 G2_zero;
    static bn128_G2 G2_one;
    static bool initialized;

    typedef bn128_Fq base_field;
    typedef bn128_Fr scalar_field;

    // Cofactor
    static const mp_size_t h_bitcount = 256;
    static const mp_size_t h_limbs = (h_bitcount+GMP_NUMB_BITS-1)/GMP_NUMB_BITS;
    static bigint<h_limbs> h;

    bn::Fp2 X, Y, Z;
    void fill_coord(bn::Fp2 coord[3]) const { coord[0] = this->X; coord[1] = this->Y; coord[2] = this->Z; };

    bn128_G2();
    bn128_G2(bn::Fp2 coord[3]) : X(coord[0]), Y(coord[1]), Z(coord[2]) {};

    void print() const;
    void print_coordinates() const;

    void to_affine_coordinates();
    void to_special();
    bool is_special() const;

    bool is_zero() const;

    bool operator==(const bn128_G2 &other) const;
    bool operator!=(const bn128_G2 &other) const;

    bn128_G2 operator+(const bn128_G2 &other) const;
    bn128_G2 operator-() const;
    bn128_G2 operator-(const bn128_G2 &other) const;

    bn128_G2 add(const bn128_G2 &other) const;
    bn128_G2 mixed_add(const bn128_G2 &other) const;
    bn128_G2 dbl() const;
    bn128_G2 mul_by_cofactor() const;

    bool is_well_formed() const;

    static bn128_G2 zero();
    static bn128_G2 one();
    static bn128_G2 random_element();

    static std::size_t size_in_bits() { return 2*base_field::ceil_size_in_bits() + 1; }
    static bigint<base_field::num_limbs> field_char() { return base_field::field_char(); }
    static bigint<scalar_field::num_limbs> order() { return scalar_field::field_char(); }

    friend std::ostream& operator<<(std::ostream &out, const bn128_G2 &g);
    friend std::istream& operator>>(std::istream &in, bn128_G2 &g);

    static void batch_to_special_all_non_zeros(std::vector<bn128_G2> &vec);
};

template<mp_size_t m>
bn128_G2 operator*(const bigint<m> &lhs, const bn128_G2 &rhs)
{
    return scalar_mul<bn128_G2, m>(rhs, lhs);
}

template<mp_size_t m, const bigint<m>& modulus_p>
bn128_G2 operator*(const Fp_model<m, modulus_p> &lhs, const bn128_G2 &rhs)
{
    return scalar_mul<bn128_G2, m>(rhs, lhs.as_bigint());
}

} // namespace libff
#endif // BN128_G2_HPP_
