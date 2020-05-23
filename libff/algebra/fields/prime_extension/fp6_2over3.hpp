/** @file
 *****************************************************************************
 Declaration of arithmetic in the finite field F[(p^3)^2]
 *****************************************************************************
 * @author     This file is part of libff, developed by SCIPR Lab
 *             and contributors (see AUTHORS).
 * @copyright  MIT license (see LICENSE file)
 *****************************************************************************/

#ifndef FP6_2OVER3_HPP_
#define FP6_2OVER3_HPP_
#include <libff/algebra/fields/prime_base/fp.hpp>
#include <libff/algebra/fields/prime_extension/fp2.hpp>
#include <libff/algebra/fields/prime_extension/fp3.hpp>

namespace libff {

/**
 * Arithmetic in the finite field F[(p^3)^2].
 *
 * Let p := modulus. This interface provides arithmetic for the extension field
 * Fp6 = Fp3[Y]/(Y^2-X) where Fp3 = Fp[X]/(X^3-non_residue) and non_residue is in Fp.
 *
 * ASSUMPTION: p = 1 (mod 6)
 */
template<mp_size_t n, const bigint<n>& modulus>
class Fp6_2over3_model;

template<mp_size_t n, const bigint<n>& modulus>
std::ostream& operator<<(std::ostream &, const Fp6_2over3_model<n, modulus> &);

template<mp_size_t n, const bigint<n>& modulus>
std::istream& operator>>(std::istream &, Fp6_2over3_model<n, modulus> &);

template<mp_size_t n, const bigint<n>& modulus>
class Fp6_2over3_model {
public:
    typedef Fp_model<n, modulus> my_Fp;
    typedef Fp2_model<n, modulus> my_Fp2;
    typedef Fp3_model<n, modulus> my_Fp3;
    typedef my_Fp3 my_Fpe;
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
    static Fp6_2over3_model<n, modulus> nqr; // a quadratic nonresidue in Fp6
    static Fp6_2over3_model<n, modulus> nqr_to_t; // nqr^t
    static my_Fp non_residue;
    static my_Fp Frobenius_coeffs_c1[6]; // non_residue^((modulus^i-1)/6)   for i=0,1,2,3,4,5

    my_Fp3 c0, c1;
    Fp6_2over3_model() {};
    Fp6_2over3_model(const my_Fp3& c0, const my_Fp3& c1) : c0(c0), c1(c1) {};

    void print() const { printf("c0/c1:\n"); c0.print(); c1.print(); }
    void clear() { c0.clear(); c1.clear(); }
    void randomize();

    bool is_zero() const { return c0.is_zero() && c1.is_zero(); }
    bool operator==(const Fp6_2over3_model &other) const;
    bool operator!=(const Fp6_2over3_model &other) const;

    Fp6_2over3_model& operator+=(const Fp6_2over3_model& other);
    Fp6_2over3_model& operator-=(const Fp6_2over3_model& other);
    Fp6_2over3_model& operator*=(const Fp6_2over3_model& other);
    Fp6_2over3_model& operator^=(const unsigned long pow);
    template<mp_size_t m>
    Fp6_2over3_model& operator^=(const bigint<m> &pow);

    Fp6_2over3_model operator+(const Fp6_2over3_model &other) const;
    Fp6_2over3_model operator-(const Fp6_2over3_model &other) const;
    Fp6_2over3_model operator*(const Fp6_2over3_model &other) const;
    Fp6_2over3_model mul_by_2345(const Fp6_2over3_model &other) const;
    Fp6_2over3_model operator^(const unsigned long pow) const;
    template<mp_size_t m>
    Fp6_2over3_model operator^(const bigint<m> &exponent) const;
    template<mp_size_t m, const bigint<m>& exp_modulus>
    Fp6_2over3_model operator^(const Fp_model<m, exp_modulus> &exponent) const;
    Fp6_2over3_model operator-() const;

    Fp6_2over3_model& square();
    Fp6_2over3_model squared() const;
    Fp6_2over3_model& invert();
    Fp6_2over3_model inverse() const;
    Fp6_2over3_model Frobenius_map(unsigned long power) const;
    Fp6_2over3_model unitary_inverse() const;
    Fp6_2over3_model cyclotomic_squared() const;
    Fp6_2over3_model sqrt() const; // HAS TO BE A SQUARE (else does not terminate)

    static my_Fp3 mul_by_non_residue(const my_Fp3 &elem);

    template<mp_size_t m>
    Fp6_2over3_model cyclotomic_exp(const bigint<m> &exponent) const;

    /** Initializes euler, s, t, t_minus_1_over_2, nqr, and nqr_to_t.
     *  Must be called before sqrt(). Alternatively, these constants can be set manually. */
    static void init_tonelli_shanks_constants();
    static std::size_t size_in_bits() { return 2*my_Fp3::size_in_bits(); }
    static constexpr std::size_t extension_degree() { return 6; }
    static constexpr bigint<n> field_char() { return modulus; }

    static Fp6_2over3_model<n, modulus> zero();
    static Fp6_2over3_model<n, modulus> one();
    static Fp6_2over3_model<n, modulus> random_element();

    friend std::ostream& operator<< <n, modulus>(std::ostream &out, const Fp6_2over3_model<n, modulus> &el);
    friend std::istream& operator>> <n, modulus>(std::istream &in, Fp6_2over3_model<n, modulus> &el);
};

#ifdef PROFILE_OP_COUNTS
template<mp_size_t n, const bigint<n>& modulus>
long long Fp6_2over3_model<n, modulus>::add_cnt = 0;

template<mp_size_t n, const bigint<n>& modulus>
long long Fp6_2over3_model<n, modulus>::sub_cnt = 0;

template<mp_size_t n, const bigint<n>& modulus>
long long Fp6_2over3_model<n, modulus>::mul_cnt = 0;

template<mp_size_t n, const bigint<n>& modulus>
long long Fp6_2over3_model<n, modulus>::sqr_cnt = 0;

template<mp_size_t n, const bigint<n>& modulus>
long long Fp6_2over3_model<n, modulus>::inv_cnt = 0;
#endif

template<mp_size_t n, const bigint<n>& modulus>
std::ostream& operator<<(std::ostream& out, const std::vector<Fp6_2over3_model<n, modulus> > &v);

template<mp_size_t n, const bigint<n>& modulus>
std::istream& operator>>(std::istream& in, std::vector<Fp6_2over3_model<n, modulus> > &v);

template<mp_size_t n, const bigint<n>& modulus>
Fp6_2over3_model<n, modulus> operator*(const Fp_model<n, modulus> &lhs, const Fp6_2over3_model<n, modulus> &rhs);

template<mp_size_t n, const bigint<n>& modulus>
bigint<6*n> Fp6_2over3_model<n, modulus>::euler;

template<mp_size_t n, const bigint<n>& modulus>
size_t Fp6_2over3_model<n, modulus>::s;

template<mp_size_t n, const bigint<n>& modulus>
bigint<6*n> Fp6_2over3_model<n, modulus>::t;

template<mp_size_t n, const bigint<n>& modulus>
bigint<6*n> Fp6_2over3_model<n, modulus>::t_minus_1_over_2;

template<mp_size_t n, const bigint<n>& modulus>
Fp6_2over3_model<n, modulus> Fp6_2over3_model<n, modulus>::nqr;

template<mp_size_t n, const bigint<n>& modulus>
Fp6_2over3_model<n, modulus> Fp6_2over3_model<n, modulus>::nqr_to_t;

template<mp_size_t n, const bigint<n>& modulus>
Fp_model<n, modulus> Fp6_2over3_model<n, modulus>::non_residue;

template<mp_size_t n, const bigint<n>& modulus>
Fp_model<n, modulus> Fp6_2over3_model<n, modulus>::Frobenius_coeffs_c1[6];

} // libff
#include <libff/algebra/fields/prime_extension/fp6_2over3.tcc>

#endif // FP6_2OVER3_HPP_
