/** @file
 *****************************************************************************

 Declaration of interfaces for the MNT4 G1 group.

 *****************************************************************************
 * @author     This file is part of libff, developed by SCIPR Lab
 *             and contributors (see AUTHORS).
 * @copyright  MIT license (see LICENSE file)
 *****************************************************************************/

#ifndef MNT4_G1_HPP_
#define MNT4_G1_HPP_

#include <vector>

#include <libff/algebra/curves/curve_utils.hpp>
#include <libff/algebra/curves/mnt/mnt4/mnt4_init.hpp>

namespace libff {

class mnt4_G1;
std::ostream& operator<<(std::ostream &, const mnt4_G1&);
std::istream& operator>>(std::istream &, mnt4_G1&);

class mnt4_G1 {
public:
#ifdef PROFILE_OP_COUNTS
    static long long add_cnt;
    static long long dbl_cnt;
#endif
    static std::vector<size_t> wnaf_window_table;
    static std::vector<size_t> fixed_base_exp_window_table;
    static mnt4_G1 G1_zero;
    static mnt4_G1 G1_one;
    static bool initialized;
    static mnt4_Fq coeff_a;
    static mnt4_Fq coeff_b;

    typedef mnt4_Fq base_field;
    typedef mnt4_Fr scalar_field;

    // Cofactor
    static const mp_size_t h_bitcount = 1;
    static const mp_size_t h_limbs = (h_bitcount+GMP_NUMB_BITS-1)/GMP_NUMB_BITS;
    static bigint<h_limbs> h;

    mnt4_Fq X, Y, Z;

    // using projective coordinates
    mnt4_G1();
    mnt4_G1(const mnt4_Fq& X, const mnt4_Fq& Y) : X(X), Y(Y), Z(base_field::one()) {}
    mnt4_G1(const mnt4_Fq& X, const mnt4_Fq& Y, const mnt4_Fq& Z) : X(X), Y(Y), Z(Z) {}

    void print() const;
    void print_coordinates() const;

    void to_affine_coordinates();
    void to_special();
    bool is_special() const;

    bool is_zero() const;

    bool operator==(const mnt4_G1 &other) const;
    bool operator!=(const mnt4_G1 &other) const;

    mnt4_G1 operator+(const mnt4_G1 &other) const;
    mnt4_G1 operator-() const;
    mnt4_G1 operator-(const mnt4_G1 &other) const;

    mnt4_G1 add(const mnt4_G1 &other) const;
    mnt4_G1 mixed_add(const mnt4_G1 &other) const;
    mnt4_G1 dbl() const;
    mnt4_G1 mul_by_cofactor() const;

    bool is_well_formed() const;

    static mnt4_G1 zero();
    static mnt4_G1 one();
    static mnt4_G1 random_element();

    static size_t size_in_bits() { return mnt4_Fq::size_in_bits() + 1; }
    static bigint<mnt4_Fq::num_limbs> base_field_char() { return mnt4_Fq::field_char(); }
    static bigint<mnt4_Fr::num_limbs> order() { return mnt4_Fr::field_char(); }

    friend std::ostream& operator<<(std::ostream &out, const mnt4_G1 &g);
    friend std::istream& operator>>(std::istream &in, mnt4_G1 &g);

    static void batch_to_special_all_non_zeros(std::vector<mnt4_G1> &vec);
};

template<mp_size_t m>
mnt4_G1 operator*(const bigint<m> &lhs, const mnt4_G1 &rhs)
{
    return scalar_mul<mnt4_G1, m>(rhs, lhs);
}

template<mp_size_t m, const bigint<m>& modulus_p>
mnt4_G1 operator*(const Fp_model<m,modulus_p> &lhs, const mnt4_G1 &rhs)
{
    return scalar_mul<mnt4_G1, m>(rhs, lhs.as_bigint());
}

std::ostream& operator<<(std::ostream& out, const std::vector<mnt4_G1> &v);
std::istream& operator>>(std::istream& in, std::vector<mnt4_G1> &v);

} // libff

#endif // MNT4_G1_HPP_
