/** @file
 *****************************************************************************

 Declaration of interfaces for the (extension) field Fp4.

 The field Fp4 equals Fp2[V]/(V^2-U) where Fp2 = Fp[U]/(U^2-non_residue) and non_residue is in Fp.

 ASSUMPTION: the modulus p is 1 mod 6.

 *****************************************************************************
 * @author     This file is part of libff, developed by SCIPR Lab
 *             and contributors (see AUTHORS).
 * @copyright  MIT license (see LICENSE file)
 *****************************************************************************/

#ifndef FP4_HPP_
#define FP4_HPP_

#include <libff/algebra/fields/prime_base/fp.hpp>
#include <libff/algebra/fields/prime_extension/fp2.hpp>

namespace libff {

template<mp_size_t n, const bigint<n>& modulus>
class Fp4_model;

template<mp_size_t n, const bigint<n>& modulus>
std::ostream& operator<<(std::ostream &, const Fp4_model<n, modulus> &);

template<mp_size_t n, const bigint<n>& modulus>
std::istream& operator>>(std::istream &, Fp4_model<n, modulus> &);

template<mp_size_t n, const bigint<n>& modulus>
class Fp4_model {
public:
    typedef Fp_model<n, modulus> my_Fp;
    typedef Fp2_model<n, modulus> my_Fp2;
    typedef my_Fp2 my_Fpe;
#ifdef PROFILE_OP_COUNTS // NOTE: op counts are affected when you exponentiate with ^
    static long long add_cnt;
    static long long sub_cnt;
    static long long mul_cnt;
    static long long sqr_cnt;
    static long long inv_cnt;
#endif

    static bigint<4*n> euler; // (modulus^4-1)/2
    static std::size_t s; // modulus^4 = 2^s * t + 1
    static bigint<4*n> t; // with t odd
    static bigint<4*n> t_minus_1_over_2; // (t-1)/2
    static Fp4_model<n, modulus> nqr; // a quadratic nonresidue in Fp4
    static Fp4_model<n, modulus> nqr_to_t; // nqr^t
    static my_Fp non_residue;
    static my_Fp Frobenius_coeffs_c1[4]; // non_residue^((modulus^i-1)/4) for i=0,1,2,3

    my_Fp2 c0, c1;
    Fp4_model() {};
    Fp4_model(const my_Fp2& c0, const my_Fp2& c1) : c0(c0), c1(c1) {};

    void print() const { printf("c0/c1:\n"); c0.print(); c1.print(); }
    void clear() { c0.clear(); c1.clear(); }
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
    bool operator==(const Fp4_model &other) const;
    bool operator!=(const Fp4_model &other) const;

    Fp4_model& operator+=(const Fp4_model& other);
    Fp4_model& operator-=(const Fp4_model& other);
    Fp4_model& operator*=(const Fp4_model& other);
    Fp4_model& operator^=(const unsigned long pow);
    template<mp_size_t m>
    Fp4_model& operator^=(const bigint<m> &pow);

    Fp4_model operator+(const Fp4_model &other) const;
    Fp4_model operator-(const Fp4_model &other) const;
    Fp4_model operator*(const Fp4_model &other) const;
    Fp4_model mul_by_023(const Fp4_model &other) const;
    Fp4_model operator^(const unsigned long pow) const;
    template<mp_size_t m>
    Fp4_model operator^(const bigint<m> &exponent) const;
    template<mp_size_t m, const bigint<m>& modulus_p>
    Fp4_model operator^(const Fp_model<m, modulus_p> &exponent) const;
    Fp4_model operator-() const;

    Fp4_model& square();
    Fp4_model squared() const;
    Fp4_model& invert();
    Fp4_model inverse() const;
    Fp4_model Frobenius_map(unsigned long power) const;
    Fp4_model unitary_inverse() const;
    Fp4_model cyclotomic_squared() const;
    Fp4_model sqrt() const; // HAS TO BE A SQUARE (else does not terminate)

    static my_Fp2 mul_by_non_residue(const my_Fp2 &elt);

    template<mp_size_t m>
    Fp4_model cyclotomic_exp(const bigint<m> &exponent) const;

    /**
     * Initializes euler, s, t, t_minus_1_over_2, nqr, and nqr_to_t.
     * Must be called before sqrt(). Alternatively, these constants can be set manually.
     */
    static void init_tonelli_shanks_constants();

    static std::size_t ceil_size_in_bits() { return 2 * my_Fp2::ceil_size_in_bits(); }
    static std::size_t floor_size_in_bits() { return 2 * my_Fp2::floor_size_in_bits(); }

    static constexpr std::size_t extension_degree() { return 4; }
    static constexpr bigint<n> field_char() { return modulus; }

    static Fp4_model<n, modulus> zero();
    static Fp4_model<n, modulus> one();
    static Fp4_model<n, modulus> random_element();

    friend std::ostream& operator<< <n, modulus>(std::ostream &out, const Fp4_model<n, modulus> &el);
    friend std::istream& operator>> <n, modulus>(std::istream &in, Fp4_model<n, modulus> &el);
};

#ifdef PROFILE_OP_COUNTS
template<mp_size_t n, const bigint<n>& modulus>
long long Fp4_model<n, modulus>::add_cnt = 0;

template<mp_size_t n, const bigint<n>& modulus>
long long Fp4_model<n, modulus>::sub_cnt = 0;

template<mp_size_t n, const bigint<n>& modulus>
long long Fp4_model<n, modulus>::mul_cnt = 0;

template<mp_size_t n, const bigint<n>& modulus>
long long Fp4_model<n, modulus>::sqr_cnt = 0;

template<mp_size_t n, const bigint<n>& modulus>
long long Fp4_model<n, modulus>::inv_cnt = 0;
#endif

template<mp_size_t n, const bigint<n>& modulus>
Fp4_model<n, modulus> operator*(const Fp_model<n, modulus> &lhs, const Fp4_model<n, modulus> &rhs);

template<mp_size_t n, const bigint<n>& modulus>
Fp4_model<n, modulus> operator*(const Fp2_model<n, modulus> &lhs, const Fp4_model<n, modulus> &rhs);

template<mp_size_t n, const bigint<n>& modulus>
bigint<4*n> Fp4_model<n, modulus>::euler;

template<mp_size_t n, const bigint<n>& modulus>
size_t Fp4_model<n, modulus>::s;

template<mp_size_t n, const bigint<n>& modulus>
bigint<4*n> Fp4_model<n, modulus>::t;

template<mp_size_t n, const bigint<n>& modulus>
bigint<4*n> Fp4_model<n, modulus>::t_minus_1_over_2;

template<mp_size_t n, const bigint<n>& modulus>
Fp4_model<n, modulus> Fp4_model<n, modulus>::nqr;

template<mp_size_t n, const bigint<n>& modulus>
Fp4_model<n, modulus> Fp4_model<n, modulus>::nqr_to_t;

template<mp_size_t n, const bigint<n>& modulus>
Fp_model<n, modulus> Fp4_model<n, modulus>::non_residue;

template<mp_size_t n, const bigint<n>& modulus>
Fp_model<n, modulus> Fp4_model<n, modulus>::Frobenius_coeffs_c1[4];


} // libff

#include <libff/algebra/fields/prime_extension/fp4.tcc>

#endif // FP4_HPP_
