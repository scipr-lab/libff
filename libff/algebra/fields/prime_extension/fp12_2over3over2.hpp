/** @file
 *****************************************************************************
 Declaration of arithmetic in the finite field F[((p^2)^3)^2].
 *****************************************************************************
 * @author     This file is part of libff, developed by SCIPR Lab
 *             and contributors (see AUTHORS).
 * @copyright  MIT license (see LICENSE file)
 *****************************************************************************/

#ifndef FP12_2OVER3OVER2_HPP_
#define FP12_2OVER3OVER2_HPP_
#include <vector>

#include <libff/algebra/fields/prime_base/fp.hpp>
#include <libff/algebra/fields/prime_extension/fp2.hpp>
#include <libff/algebra/fields/prime_extension/fp6_3over2.hpp>

namespace libff {

template<mp_size_t n, const bigint<n>& modulus>
class Fp12_2over3over2_model;

template<mp_size_t n, const bigint<n>& modulus>
std::ostream& operator<<(std::ostream &, const Fp12_2over3over2_model<n, modulus> &);

template<mp_size_t n, const bigint<n>& modulus>
std::istream& operator>>(std::istream &, Fp12_2over3over2_model<n, modulus> &);

/**
 * Arithmetic in the finite field F[((p^2)^3)^2].
 *
 * Let p := modulus. This interface provides arithmetic for the extension field
 * Fp12 = Fp6[W]/(W^2-V) where Fp6 = Fp2[V]/(V^3-non_residue) and non_residue is in Fp2
 *
 * ASSUMPTION: p = 1 (mod 6)
 */
template<mp_size_t n, const bigint<n>& modulus>
class Fp12_2over3over2_model {
public:
    typedef Fp_model<n, modulus> my_Fp;
    typedef Fp2_model<n, modulus> my_Fp2;
    typedef Fp6_3over2_model<n, modulus> my_Fp6;
#ifdef PROFILE_OP_COUNTS // NOTE: op counts are affected when you exponentiate with ^
    static long long add_cnt;
    static long long sub_cnt;
    static long long mul_cnt;
    static long long sqr_cnt;
    static long long inv_cnt;
#endif

    static bigint<12*n> euler; // (modulus^12-1)/2
    static std::size_t s; // modulus^12 = 2^s * t + 1
    static bigint<12*n> t; // with t odd
    static bigint<12*n> t_minus_1_over_2; // (t-1)/2
    static Fp12_2over3over2_model<n, modulus> nqr; // a quadratic nonresidue in Fp12
    static Fp12_2over3over2_model<n, modulus> nqr_to_t; // nqr^t
    static Fp2_model<n, modulus> non_residue;
    static Fp2_model<n, modulus> Frobenius_coeffs_c1[12]; // non_residue^((modulus^i-1)/6) for i=0,...,11

    my_Fp6 c0, c1;
    Fp12_2over3over2_model() {};
    Fp12_2over3over2_model(const my_Fp6& c0, const my_Fp6& c1) : c0(c0), c1(c1) {};

    void clear() { c0.clear(); c1.clear(); }
    void print() const { printf("c0/c1:\n"); c0.print(); c1.print(); }
    void randomize();

    /**
     * Returns the constituent bits in 64 bit words, in little-endian order.
     * Only the right-most ceil_size_in_bits() bits are used; other bits are 0.
     */
    std::vector<uint64_t> to_words() const;
    /**
     * Sets the field element from the given bits in 64 bit words, in little-endian order.
     * Only the right-most ceil_size_in_bits() bits are used; other bits are ignored.
     * Returns true when the right-most bits of each element represent a value less than the modulus.
     */
    bool from_words(std::vector<uint64_t> words);

    bool is_zero() const { return c0.is_zero() && c1.is_zero(); }
    bool operator==(const Fp12_2over3over2_model &other) const;
    bool operator!=(const Fp12_2over3over2_model &other) const;

    Fp12_2over3over2_model& operator+=(const Fp12_2over3over2_model& other);
    Fp12_2over3over2_model& operator-=(const Fp12_2over3over2_model& other);
    Fp12_2over3over2_model& operator*=(const Fp12_2over3over2_model& other);
    Fp12_2over3over2_model& operator^=(const unsigned long pow);
    template<mp_size_t m>
    Fp12_2over3over2_model& operator^=(const bigint<m> &pow);

    Fp12_2over3over2_model operator+(const Fp12_2over3over2_model &other) const;
    Fp12_2over3over2_model operator-(const Fp12_2over3over2_model &other) const;
    Fp12_2over3over2_model operator*(const Fp12_2over3over2_model &other) const;
    Fp12_2over3over2_model operator^(const unsigned long pow) const;
    template<mp_size_t m>
    Fp12_2over3over2_model operator^(const bigint<m> &exponent) const;
    template<mp_size_t m, const bigint<m>& exp_modulus>
    Fp12_2over3over2_model operator^(const Fp_model<m, exp_modulus> &exponent) const;
    Fp12_2over3over2_model operator-() const;

    Fp12_2over3over2_model& square();
    Fp12_2over3over2_model squared() const; // default is squared_complex
    Fp12_2over3over2_model squared_karatsuba() const;
    Fp12_2over3over2_model squared_complex() const;
    Fp12_2over3over2_model& invert();
    Fp12_2over3over2_model inverse() const;
    Fp12_2over3over2_model Frobenius_map(unsigned long power) const;
    Fp12_2over3over2_model unitary_inverse() const;
    Fp12_2over3over2_model cyclotomic_squared() const;
    Fp12_2over3over2_model sqrt() const; // HAS TO BE A SQUARE (else does not terminate)

    Fp12_2over3over2_model mul_by_024(const my_Fp2 &ell_0, const my_Fp2 &ell_VW, const my_Fp2 &ell_VV) const;
    Fp12_2over3over2_model mul_by_045(const my_Fp2 &ell_0, const my_Fp2 &ell_VW, const my_Fp2 &ell_VV) const;

    static my_Fp6 mul_by_non_residue(const my_Fp6 &elt);

    template<mp_size_t m>
    Fp12_2over3over2_model cyclotomic_exp(const bigint<m> &exponent) const;

    static std::size_t ceil_size_in_bits() { return 2 * my_Fp6::ceil_size_in_bits(); }
    static std::size_t floor_size_in_bits() { return 2 * my_Fp6::floor_size_in_bits(); }

    static constexpr std::size_t extension_degree() { return 12; }
    static constexpr bigint<n> field_char() { return modulus; }

    static Fp12_2over3over2_model<n, modulus> zero();
    static Fp12_2over3over2_model<n, modulus> one();
    static Fp12_2over3over2_model<n, modulus> random_element();

    friend std::ostream& operator<< <n, modulus>(std::ostream &out, const Fp12_2over3over2_model<n, modulus> &el);
    friend std::istream& operator>> <n, modulus>(std::istream &in, Fp12_2over3over2_model<n, modulus> &el);
};

#ifdef PROFILE_OP_COUNTS
template<mp_size_t n, const bigint<n>& modulus>
long long Fp12_2over3over2_model<n, modulus>::add_cnt = 0;

template<mp_size_t n, const bigint<n>& modulus>
long long Fp12_2over3over2_model<n, modulus>::sub_cnt = 0;

template<mp_size_t n, const bigint<n>& modulus>
long long Fp12_2over3over2_model<n, modulus>::mul_cnt = 0;

template<mp_size_t n, const bigint<n>& modulus>
long long Fp12_2over3over2_model<n, modulus>::sqr_cnt = 0;

template<mp_size_t n, const bigint<n>& modulus>
long long Fp12_2over3over2_model<n, modulus>::inv_cnt = 0;
#endif

template<mp_size_t n, const bigint<n>& modulus>
std::ostream& operator<<(std::ostream& out, const std::vector<Fp12_2over3over2_model<n, modulus> > &v);

template<mp_size_t n, const bigint<n>& modulus>
std::istream& operator>>(std::istream& in, std::vector<Fp12_2over3over2_model<n, modulus> > &v);

template<mp_size_t n, const bigint<n>& modulus>
Fp12_2over3over2_model<n, modulus> operator*(const Fp_model<n, modulus> &lhs, const Fp12_2over3over2_model<n, modulus> &rhs);

template<mp_size_t n, const bigint<n>& modulus>
Fp12_2over3over2_model<n, modulus> operator*(const Fp2_model<n, modulus> &lhs, const Fp12_2over3over2_model<n, modulus> &rhs);

template<mp_size_t n, const bigint<n>& modulus>
Fp12_2over3over2_model<n, modulus> operator*(const Fp6_3over2_model<n, modulus> &lhs, const Fp12_2over3over2_model<n, modulus> &rhs);

template<mp_size_t n, const bigint<n>& modulus>
bigint<12*n> Fp12_2over3over2_model<n, modulus>::euler;

template<mp_size_t n, const bigint<n>& modulus>
size_t Fp12_2over3over2_model<n, modulus>::s;

template<mp_size_t n, const bigint<n>& modulus>
bigint<12*n> Fp12_2over3over2_model<n, modulus>::t;

template<mp_size_t n, const bigint<n>& modulus>
bigint<12*n> Fp12_2over3over2_model<n, modulus>::t_minus_1_over_2;

template<mp_size_t n, const bigint<n>& modulus>
Fp12_2over3over2_model<n, modulus> Fp12_2over3over2_model<n, modulus>::nqr;

template<mp_size_t n, const bigint<n>& modulus>
Fp12_2over3over2_model<n, modulus> Fp12_2over3over2_model<n, modulus>::nqr_to_t;

template<mp_size_t n, const bigint<n>& modulus>
Fp2_model<n, modulus> Fp12_2over3over2_model<n, modulus>::non_residue;

template<mp_size_t n, const bigint<n>& modulus>
Fp2_model<n, modulus> Fp12_2over3over2_model<n, modulus>::Frobenius_coeffs_c1[12];

} // namespace libff
#include <libff/algebra/fields/prime_extension/fp12_2over3over2.tcc>
#endif // FP12_2OVER3OVER2_HPP_
