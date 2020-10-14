/** @file
 *****************************************************************************
 Declaration of arithmetic in the finite field F[(p^2)^3]
 *****************************************************************************
 * @author     This file is part of libff, developed by SCIPR Lab
 *             and contributors (see AUTHORS).
 * @copyright  MIT license (see LICENSE file)
 *****************************************************************************/

#ifndef FP6_3OVER2_HPP_
#define FP6_3OVER2_HPP_
#include <vector>

#include <libff/algebra/fields/prime_base/fp.hpp>
#include <libff/algebra/fields/prime_extension/fp2.hpp>

namespace libff {

template<mp_size_t n, const bigint<n>& modulus>
class Fp6_3over2_model;

template<mp_size_t n, const bigint<n>& modulus>
std::ostream& operator<<(std::ostream &, const Fp6_3over2_model<n, modulus> &);

template<mp_size_t n, const bigint<n>& modulus>
std::istream& operator>>(std::istream &, Fp6_3over2_model<n, modulus> &);

/**
 * Arithmetic in the finite field F[(p^2)^3].
 *
 * Let p := modulus. This interface provides arithmetic for the extension field
 *  Fp6 = Fp2[V]/(V^3-non_residue) where non_residue is in Fp.
 *
 * ASSUMPTION: p = 1 (mod 6)
 */
template<mp_size_t n, const bigint<n>& modulus>
class Fp6_3over2_model {
public:
    typedef Fp_model<n, modulus> my_Fp;
    typedef Fp2_model<n, modulus> my_Fp2;
#ifdef PROFILE_OP_COUNTS // NOTE: op counts are affected when you exponentiate with ^
    static long long add_cnt;
    static long long sub_cnt;
    static long long mul_cnt;
    static long long sqr_cnt;
    static long long inv_cnt;
#endif

    static bigint<6*n> euler; // (modulus^6-1)/2
    static std::size_t s; // modulus^6 = 2^s * t + 1
    static bigint<6*n> t; // with t odd
    static bigint<6*n> t_minus_1_over_2; // (t-1)/2
    static Fp6_3over2_model<n, modulus> nqr; // a quadratic nonresidue in Fp6
    static Fp6_3over2_model<n, modulus> nqr_to_t; // nqr^t
    static my_Fp2 non_residue;
    static my_Fp2 Frobenius_coeffs_c1[6]; // non_residue^((modulus^i-1)/3)   for i=0,1,2,3,4,5
    static my_Fp2 Frobenius_coeffs_c2[6]; // non_residue^((2*modulus^i-2)/3) for i=0,1,2,3,4,5

    my_Fp2 c0, c1, c2;
    Fp6_3over2_model() {};
    Fp6_3over2_model(const my_Fp2& c0, const my_Fp2& c1, const my_Fp2& c2) : c0(c0), c1(c1), c2(c2) {};

    void clear() { c0.clear(); c1.clear(); c2.clear(); }
    void print() const { printf("c0/c1/c2:\n"); c0.print(); c1.print(); c2.print(); }
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

    bool is_zero() const { return c0.is_zero() && c1.is_zero() && c2.is_zero(); }
    bool operator==(const Fp6_3over2_model &other) const;
    bool operator!=(const Fp6_3over2_model &other) const;

    Fp6_3over2_model& operator+=(const Fp6_3over2_model& other);
    Fp6_3over2_model& operator-=(const Fp6_3over2_model& other);
    Fp6_3over2_model& operator*=(const Fp6_3over2_model& other);
    Fp6_3over2_model& operator^=(const unsigned long pow);
    template<mp_size_t m>
    Fp6_3over2_model& operator^=(const bigint<m> &pow);

    Fp6_3over2_model operator+(const Fp6_3over2_model &other) const;
    Fp6_3over2_model operator-(const Fp6_3over2_model &other) const;
    Fp6_3over2_model operator*(const Fp6_3over2_model &other) const;
    Fp6_3over2_model operator^(const unsigned long pow) const;
    template<mp_size_t m>
    Fp6_3over2_model operator^(const bigint<m> &other) const;
    Fp6_3over2_model operator-() const;

    Fp6_3over2_model& square();
    Fp6_3over2_model squared() const;
    Fp6_3over2_model& invert();
    Fp6_3over2_model inverse() const;
    Fp6_3over2_model Frobenius_map(unsigned long power) const;
    Fp6_3over2_model sqrt() const; // HAS TO BE A SQUARE (else does not terminate)

    static my_Fp2 mul_by_non_residue(const my_Fp2 &elt);

    /**
     * Initializes euler, s, t, t_minus_1_over_2, nqr, and nqr_to_t.
     * Must be called before sqrt(). Alternatively, these constants can be set manually.
     */
    static void init_tonelli_shanks_constants();

    static std::size_t ceil_size_in_bits() { return 3 * my_Fp2::ceil_size_in_bits(); }
    static std::size_t floor_size_in_bits() { return 3 * my_Fp2::floor_size_in_bits(); }

    static constexpr std::size_t extension_degree() { return 6; }
    static constexpr bigint<n> field_char() { return modulus; }

    static Fp6_3over2_model<n, modulus> zero();
    static Fp6_3over2_model<n, modulus> one();
    static Fp6_3over2_model<n, modulus> random_element();

    friend std::ostream& operator<< <n, modulus>(std::ostream &out, const Fp6_3over2_model<n, modulus> &el);
    friend std::istream& operator>> <n, modulus>(std::istream &in, Fp6_3over2_model<n, modulus> &el);
};

#ifdef PROFILE_OP_COUNTS
template<mp_size_t n, const bigint<n>& modulus>
long long Fp6_3over2_model<n, modulus>::add_cnt = 0;

template<mp_size_t n, const bigint<n>& modulus>
long long Fp6_3over2_model<n, modulus>::sub_cnt = 0;

template<mp_size_t n, const bigint<n>& modulus>
long long Fp6_3over2_model<n, modulus>::mul_cnt = 0;

template<mp_size_t n, const bigint<n>& modulus>
long long Fp6_3over2_model<n, modulus>::sqr_cnt = 0;

template<mp_size_t n, const bigint<n>& modulus>
long long Fp6_3over2_model<n, modulus>::inv_cnt = 0;
#endif

template<mp_size_t n, const bigint<n>& modulus>
std::ostream& operator<<(std::ostream& out, const std::vector<Fp6_3over2_model<n, modulus> > &v);

template<mp_size_t n, const bigint<n>& modulus>
std::istream& operator>>(std::istream& in, std::vector<Fp6_3over2_model<n, modulus> > &v);

template<mp_size_t n, const bigint<n>& modulus>
Fp6_3over2_model<n, modulus> operator*(const Fp_model<n, modulus> &lhs, const Fp6_3over2_model<n, modulus> &rhs);

template<mp_size_t n, const bigint<n>& modulus>
Fp6_3over2_model<n, modulus> operator*(const Fp2_model<n, modulus> &lhs, const Fp6_3over2_model<n, modulus> &rhs);

template<mp_size_t n, const bigint<n>& modulus>
bigint<6*n> Fp6_3over2_model<n, modulus>::euler;

template<mp_size_t n, const bigint<n>& modulus>
size_t Fp6_3over2_model<n, modulus>::s;

template<mp_size_t n, const bigint<n>& modulus>
bigint<6*n> Fp6_3over2_model<n, modulus>::t;

template<mp_size_t n, const bigint<n>& modulus>
bigint<6*n> Fp6_3over2_model<n, modulus>::t_minus_1_over_2;

template<mp_size_t n, const bigint<n>& modulus>
Fp6_3over2_model<n, modulus> Fp6_3over2_model<n, modulus>::nqr;

template<mp_size_t n, const bigint<n>& modulus>
Fp6_3over2_model<n, modulus> Fp6_3over2_model<n, modulus>::nqr_to_t;

template<mp_size_t n, const bigint<n>& modulus>
Fp2_model<n, modulus> Fp6_3over2_model<n, modulus>::non_residue;

template<mp_size_t n, const bigint<n>& modulus>
Fp2_model<n, modulus> Fp6_3over2_model<n, modulus>::Frobenius_coeffs_c1[6];

template<mp_size_t n, const bigint<n>& modulus>
Fp2_model<n, modulus> Fp6_3over2_model<n, modulus>::Frobenius_coeffs_c2[6];

} // libff
#include <libff/algebra/fields/prime_extension/fp6_3over2.tcc>

#endif // FP6_3OVER2_HPP_
