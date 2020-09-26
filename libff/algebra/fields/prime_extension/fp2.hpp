/** @file
 *****************************************************************************
 Implementation of arithmetic in the finite field F[p^2].
 *****************************************************************************
 * @author     This file is part of libff, developed by SCIPR Lab
 *             and contributors (see AUTHORS).
 * @copyright  MIT license (see LICENSE file)
 *****************************************************************************/

#ifndef FP2_HPP_
#define FP2_HPP_
#include <vector>

#include <libff/algebra/fields/prime_base/fp.hpp>

namespace libff {

template<mp_size_t n, const bigint<n>& modulus>
class Fp2_model;

template<mp_size_t n, const bigint<n>& modulus>
std::ostream& operator<<(std::ostream &, const Fp2_model<n, modulus> &);

template<mp_size_t n, const bigint<n>& modulus>
std::istream& operator>>(std::istream &, Fp2_model<n, modulus> &);

/**
 * Arithmetic in the field F[p^2].
 *
 * Let p := modulus. This interface provides arithmetic for the extension field
 * Fp2 = Fp[U]/(U^2-non_residue), where non_residue is in Fp.
 *
 * ASSUMPTION: p = 1 (mod 6)
 */
template<mp_size_t n, const bigint<n>& modulus>
class Fp2_model {
public:
    typedef Fp_model<n, modulus> my_Fp;
#ifdef PROFILE_OP_COUNTS // NOTE: op counts are affected when you exponentiate with ^
    static long long add_cnt;
    static long long sub_cnt;
    static long long mul_cnt;
    static long long sqr_cnt;
    static long long inv_cnt;
#endif

    static bigint<2*n> euler; // (modulus^2-1)/2
    static std::size_t s;       // modulus^2 = 2^s * t + 1
    static bigint<2*n> t;  // with t odd
    static bigint<2*n> t_minus_1_over_2; // (t-1)/2
    static my_Fp non_residue; // X^4-non_residue irreducible over Fp; used for constructing Fp2 = Fp[X] / (X^2 - non_residue)
    static Fp2_model<n, modulus> nqr; // a quadratic nonresidue in Fp2
    static Fp2_model<n, modulus> nqr_to_t; // nqr^t
    static my_Fp Frobenius_coeffs_c1[2]; // non_residue^((modulus^i-1)/2) for i=0,1

    my_Fp c0, c1;
    Fp2_model() {};
    Fp2_model(const my_Fp& c0, const my_Fp& c1) : c0(c0), c1(c1) {};

    void clear() { c0.clear(); c1.clear(); }
    void print() const { printf("c0/c1:\n"); c0.print(); c1.print(); }
    void randomize();

    bool is_zero() const { return c0.is_zero() && c1.is_zero(); }
    bool operator==(const Fp2_model &other) const;
    bool operator!=(const Fp2_model &other) const;

    Fp2_model& operator+=(const Fp2_model& other);
    Fp2_model& operator-=(const Fp2_model& other);
    Fp2_model& operator*=(const Fp2_model& other);
    Fp2_model& operator^=(const unsigned long pow);
    template<mp_size_t m>
    Fp2_model& operator^=(const bigint<m> &pow);

    Fp2_model operator+(const Fp2_model &other) const;
    Fp2_model operator-(const Fp2_model &other) const;
    Fp2_model operator*(const Fp2_model &other) const;
    Fp2_model operator^(const unsigned long pow) const;
    template<mp_size_t m>
    Fp2_model operator^(const bigint<m> &other) const;
    Fp2_model operator-() const;

    Fp2_model& square(); // default is squared_complex
    Fp2_model squared() const; // default is squared_complex
    Fp2_model& invert();
    Fp2_model inverse() const;
    Fp2_model Frobenius_map(unsigned long power) const;
    Fp2_model sqrt() const; // HAS TO BE A SQUARE (else does not terminate)
    Fp2_model squared_karatsuba() const;
    Fp2_model squared_complex() const;

    /** Initializes euler, s, t, t_minus_1_over_2, nqr, and nqr_to_t.
     *  Must be called before sqrt(). Alternatively, these constants can be set manually. */
    static void init_tonelli_shanks_constants();
    static std::size_t ceil_size_in_bits() { return 2*my_Fp::ceil_size_in_bits(); }
    static constexpr std::size_t extension_degree() { return 2; }
    static constexpr bigint<n> field_char() { return modulus; }

    static Fp2_model<n, modulus> zero();
    static Fp2_model<n, modulus> one();
    static Fp2_model<n, modulus> random_element();

    friend std::ostream& operator<< <n, modulus>(std::ostream &out, const Fp2_model<n, modulus> &el);
    friend std::istream& operator>> <n, modulus>(std::istream &in, Fp2_model<n, modulus> &el);
};

#ifdef PROFILE_OP_COUNTS
template<mp_size_t n, const bigint<n>& modulus>
long long Fp2_model<n, modulus>::add_cnt = 0;

template<mp_size_t n, const bigint<n>& modulus>
long long Fp2_model<n, modulus>::sub_cnt = 0;

template<mp_size_t n, const bigint<n>& modulus>
long long Fp2_model<n, modulus>::mul_cnt = 0;

template<mp_size_t n, const bigint<n>& modulus>
long long Fp2_model<n, modulus>::sqr_cnt = 0;

template<mp_size_t n, const bigint<n>& modulus>
long long Fp2_model<n, modulus>::inv_cnt = 0;
#endif

template<mp_size_t n, const bigint<n>& modulus>
std::ostream& operator<<(std::ostream& out, const std::vector<Fp2_model<n, modulus> > &v);

template<mp_size_t n, const bigint<n>& modulus>
std::istream& operator>>(std::istream& in, std::vector<Fp2_model<n, modulus> > &v);

template<mp_size_t n, const bigint<n>& modulus>
Fp2_model<n, modulus> operator*(const Fp_model<n, modulus> &lhs, const Fp2_model<n, modulus> &rhs);

template<mp_size_t n, const bigint<n>& modulus>
bigint<2*n> Fp2_model<n, modulus>::euler;

template<mp_size_t n, const bigint<n>& modulus>
size_t Fp2_model<n, modulus>::s;

template<mp_size_t n, const bigint<n>& modulus>
bigint<2*n> Fp2_model<n, modulus>::t;

template<mp_size_t n, const bigint<n>& modulus>
bigint<2*n> Fp2_model<n, modulus>::t_minus_1_over_2;

template<mp_size_t n, const bigint<n>& modulus>
Fp_model<n, modulus> Fp2_model<n, modulus>::non_residue;

template<mp_size_t n, const bigint<n>& modulus>
Fp2_model<n, modulus> Fp2_model<n, modulus>::nqr;

template<mp_size_t n, const bigint<n>& modulus>
Fp2_model<n, modulus> Fp2_model<n, modulus>::nqr_to_t;

template<mp_size_t n, const bigint<n>& modulus>
Fp_model<n, modulus> Fp2_model<n, modulus>::Frobenius_coeffs_c1[2];

} // libff
#include <libff/algebra/fields/prime_extension/fp2.tcc>

#endif // FP2_HPP_
