/** @file
 *****************************************************************************

 Declaration of interfaces for the MNT4 G2 group.

 *****************************************************************************
 * @author     This file is part of libff, developed by SCIPR Lab
 *             and contributors (see AUTHORS).
 * @copyright  MIT license (see LICENSE file)
 *****************************************************************************/

#ifndef MNT4_G2_HPP_
#define MNT4_G2_HPP_

#include <vector>

#include <libff/algebra/curves/curve_utils.hpp>
#include <libff/algebra/curves/mnt/mnt4/mnt4_init.hpp>

namespace libff {

class mnt4_G2;
std::ostream& operator<<(std::ostream &, const mnt4_G2&);
std::istream& operator>>(std::istream &, mnt4_G2&);

class mnt4_G2 {
public:
#ifdef PROFILE_OP_COUNTS
    static long long add_cnt;
    static long long dbl_cnt;
#endif
    static std::vector<size_t> wnaf_window_table;
    static std::vector<size_t> fixed_base_exp_window_table;
    static mnt4_G2 G2_zero;
    static mnt4_G2 G2_one;
    static bool initialized;
    static mnt4_Fq2 twist;
    static mnt4_Fq2 coeff_a;
    static mnt4_Fq2 coeff_b;

    typedef mnt4_Fq base_field;
    typedef mnt4_Fq2 twist_field;
    typedef mnt4_Fr scalar_field;

    // Cofactor
    static const mp_size_t h_bitcount = 298;
    static const mp_size_t h_limbs = (h_bitcount+GMP_NUMB_BITS-1)/GMP_NUMB_BITS;
    static bigint<h_limbs> h;

    mnt4_Fq2 X, Y, Z;

    // using projective coordinates
    mnt4_G2();
    mnt4_G2(const mnt4_Fq2& X, const mnt4_Fq2& Y, const mnt4_Fq2& Z) : X(X), Y(Y), Z(Z) {};

    static mnt4_Fq2 mul_by_a(const mnt4_Fq2 &elt);
    static mnt4_Fq2 mul_by_b(const mnt4_Fq2 &elt);

    void print() const;
    void print_coordinates() const;

    void to_affine_coordinates();
    void to_special();
    bool is_special() const;

    bool is_zero() const;

    bool operator==(const mnt4_G2 &other) const;
    bool operator!=(const mnt4_G2 &other) const;

    mnt4_G2 operator+(const mnt4_G2 &other) const;
    mnt4_G2 operator-() const;
    mnt4_G2 operator-(const mnt4_G2 &other) const;

    mnt4_G2 add(const mnt4_G2 &other) const;
    mnt4_G2 mixed_add(const mnt4_G2 &other) const;
    mnt4_G2 dbl() const;
    mnt4_G2 mul_by_q() const;
    mnt4_G2 mul_by_cofactor() const;

    bool is_well_formed() const;

    static mnt4_G2 zero();
    static mnt4_G2 one();
    static mnt4_G2 random_element();

    static size_t size_in_bits() { return mnt4_Fq2::size_in_bits() + 1; }
    static bigint<mnt4_Fq::num_limbs> base_field_char() { return mnt4_Fq::field_char(); }
    static bigint<mnt4_Fr::num_limbs> order() { return mnt4_Fr::field_char(); }

    friend std::ostream& operator<<(std::ostream &out, const mnt4_G2 &g);
    friend std::istream& operator>>(std::istream &in, mnt4_G2 &g);

    static void batch_to_special_all_non_zeros(std::vector<mnt4_G2> &vec);
};

template<mp_size_t m>
mnt4_G2 operator*(const bigint<m> &lhs, const mnt4_G2 &rhs)
{
    return scalar_mul<mnt4_G2, m>(rhs, lhs);
}

template<mp_size_t m, const bigint<m>& modulus_p>
mnt4_G2 operator*(const Fp_model<m,modulus_p> &lhs, const mnt4_G2 &rhs)
{
    return scalar_mul<mnt4_G2, m>(rhs, lhs.as_bigint());
}

} // libff

#endif // MNT4_G2_HPP_
