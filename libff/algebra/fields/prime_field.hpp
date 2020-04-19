/** @file
 *****************************************************************************
 Declaration of common API for all finite fields in the prime_base/ and
 prime_extension/ directories.

 Includes fields Fp^n for specified n. All of the prime extension fields must
 implement all functions declared in this class.

 However, this class is not actually the parent class of any field. All APIs
 are enforced through tests instead.
 *****************************************************************************
 * @author     This file is part of libff, developed by SCIPR Lab
 *             and contributors (see AUTHORS).
 * @copyright  MIT license (see LICENSE file)
 *****************************************************************************/
#include <libff/algebra/fields/bigint.hpp>
#include <vector>

namespace libff {

template<typename T, mp_size_t n, const bigint<n>& modulus>
class PrimeField;

/* The type parameter T is intended to be set to the child class
 * when this class is extended. For example,
 * class Fp_model : public PrimeField<Fp_model, n, modulus> ...
 */
template<typename T, mp_size_t n, const bigint<n>& modulus>
class PrimeField {
    /* Functions unique to prime fields */

    /** If extension field, returns the base field's characteristic. */
    static constexpr bigint<n> field_char(); // has not been implemented in Fp or gf2^n

    /** If base field, is the identity. */
    T Frobenius_map(unsigned long power) const;

    /* Functions common to all finite fields */

    virtual T& operator+=(const T& other) = 0; // has not been implemented in fp2 and above
    virtual T& operator-=(const T& other) = 0; // has not been implemented in fp2 and above
    virtual T& operator*=(const T& other) = 0; // has not been implemented in fp2 and above
    virtual T& operator^=(const unsigned long pow) = 0; // has not been implemented in gf2^n or fp2 and above
    template<mp_size_t m>
    virtual T& operator^=(const bigint<m> &pow) = 0; // has not been implemented in gf2^n or fp2 and above

    virtual T& square() = 0; // has not been implemented in Fp^n
    virtual T& invert() = 0; // has not been implemented in gf2^n or fp2 and above

    virtual T operator+(const T& other) const;
    virtual T operator-(const T& other) const;
    virtual T operator*(const T& other) const;
    virtual T operator^(const unsigned long pow) const; // has not been implemented in gf2^n
    template<mp_size_t m>
    virtual T operator^(const bigint<m> &pow) const; // has not been implemented in gf2^n
    virtual T operator-() const = 0;

    virtual T squared() const;
    virtual T inverse() const;
    /** HAS TO BE A SQUARE (else does not terminate). */
    virtual T sqrt() const = 0; // has not been implemented in gf2^n or fp4 and above

    bool operator==(const T& other) const = 0;
    bool operator!=(const T& other) const = 0;
    bool is_zero() const = 0;

    void print() const = 0;

    void randomize() = 0; // has not been implemented in Fp^n
    void clear() = 0; // has not been implemented in gf2^n

    // the following should be defined in child classes, but are static so they can't be inherited
    static T zero();
    static T one();
    static T random_element();
    /** Equals 1 for prime field Fp. */
    static constexpr std::size_t extension_degree();
    static std::size_t size_in_bits();

    // the following should be defined as well but can't be inherited
    friend std::ostream& operator<< <n,modulus>(std::ostream &out, const Fp_model<n, modulus> &p); // has not been implemented in gf2^n
    friend std::istream& operator>> <n,modulus>(std::istream &in, Fp_model<n, modulus> &p); // has not been implemented in gf2^n
};

// has not been implemented in gf2^n
#ifdef PROFILE_OP_COUNTS
template<typename T>
long long Field<T>::add_cnt = 0;

template<typename T>
long long Field<T>::sub_cnt = 0;

template<typename T>
long long Field<T>::mul_cnt = 0;

template<typename T>
long long Field<T>::sqr_cnt = 0;

template<typename T>
long long Field<T>::inv_cnt = 0;
#endif

} // libff
