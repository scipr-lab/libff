/** @file
 *****************************************************************************
 * @author     This file is part of libff, developed by SCIPR Lab
 *             and contributors (see AUTHORS).
 * @copyright  MIT license (see LICENSE file)
 *****************************************************************************/

#ifndef EDWARDS_G1_HPP_
#define EDWARDS_G1_HPP_
#include <vector>

#include <libff/algebra/curves/curve_utils.hpp>
#include <libff/algebra/curves/edwards/edwards_init.hpp>

namespace libff {

class edwards_G1;
std::ostream& operator<<(std::ostream &, const edwards_G1&);
std::istream& operator>>(std::istream &, edwards_G1&);

class edwards_G1 {
public:
#ifdef PROFILE_OP_COUNTS
    static long long add_cnt;
    static long long dbl_cnt;
#endif
    static std::vector<size_t> wnaf_window_table;
    static std::vector<size_t> fixed_base_exp_window_table;
    static edwards_G1 G1_zero;
    static edwards_G1 G1_one;
    static bool initialized;

    edwards_Fq X, Y, Z;
    edwards_G1();
private:
    edwards_G1(const edwards_Fq& X, const edwards_Fq& Y, const edwards_Fq& Z) : X(X), Y(Y), Z(Z) {};

public:
    typedef edwards_Fq base_field;
    typedef edwards_Fr scalar_field;
    // using inverted coordinates
    edwards_G1(const edwards_Fq& X, const edwards_Fq& Y) : X(Y), Y(X), Z(X*Y) {};

    void print() const;
    void print_coordinates() const;

    void to_affine_coordinates();
    void to_special();
    bool is_special() const;

    bool is_zero() const;

    bool operator==(const edwards_G1 &other) const;
    bool operator!=(const edwards_G1 &other) const;

    edwards_G1 operator+(const edwards_G1 &other) const;
    edwards_G1 operator-() const;
    edwards_G1 operator-(const edwards_G1 &other) const;

    edwards_G1 add(const edwards_G1 &other) const;
    edwards_G1 mixed_add(const edwards_G1 &other) const;
    edwards_G1 dbl() const;

    bool is_well_formed() const;

    static edwards_G1 zero();
    static edwards_G1 one();
    static edwards_G1 random_element();

    static size_t size_in_bits() { return edwards_Fq::size_in_bits() + 1; }
    static bigint<base_field::num_limbs> base_field_char() { return base_field::field_char(); }
    static bigint<scalar_field::num_limbs> order() { return scalar_field::field_char(); }

    void write_uncompressed(std::ostream &) const;
    void write_compressed(std::ostream &) const;
    static void read_uncompressed(std::istream &, edwards_G1 &);
    static void read_compressed(std::istream &, edwards_G1 &);

    static void batch_to_special_all_non_zeros(std::vector<edwards_G1> &vec);
};

template<mp_size_t m>
edwards_G1 operator*(const bigint<m> &lhs, const edwards_G1 &rhs)
{
    return scalar_mul<edwards_G1, m>(rhs, lhs);
}

template<mp_size_t m, const bigint<m>& modulus_p>
edwards_G1 operator*(const Fp_model<m,modulus_p> &lhs, const edwards_G1 &rhs)
{
    return scalar_mul<edwards_G1, m>(rhs, lhs.as_bigint());
}

} // libff
#endif // EDWARDS_G1_HPP_
