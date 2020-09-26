/** @file
 *****************************************************************************
 Declaration of common API for all finite fields in the prime_base/ and
 prime_extension/ directories.

 Currently NOT used by the fields in this library. This class is not actually
 the parent class of any field. All APIs are enforced through tests instead.

 The reason for this is to ensure high performance of all fields. This class
 exists as documentation for common API between fields.

 Includes fields Fp^n for specified n. All of the prime extension fields must
 implement all functions declared in this class.
 *****************************************************************************
 * @author     This file is part of libff, developed by SCIPR Lab
 *             and contributors (see AUTHORS).
 * @copyright  MIT license (see LICENSE file)
 *****************************************************************************/
#include <libff/algebra/field_utils/bigint.hpp>
#include <vector>

namespace libff {

template<typename T, mp_size_t n, const bigint<n>& modulus>
class FpnField;

/* The type parameter T is intended to be set to the child class
 * when this class is extended. For example,
 * class Fp_model : public FpnField<Fp_model, n, modulus> ...
 */
template<typename T, mp_size_t n, const bigint<n>& modulus>
class FpnField {
public:
    /* Functions unique to Fp^n fields */

    /** The following parameters are used for the Tonelli-Shanks square root algorithms.
     *  Can be set manually or calculated with init_tonelli_shanks_constants(). */
    static bigint<T::extension_degree()*n> euler; // (modulus^deg-1)/2
    static std::size_t s; // modulus^deg = 2^s * t + 1
    static bigint<T::extension_degree()*n> t; // with t odd
    static bigint<T::extension_degree()*n> t_minus_1_over_2; // (t-1)/2
    static T nqr; // a quadratic nonresidue in field
    static T nqr_to_t; // nqr^t

    /** If extension field, returns the base field's characteristic. */
    static constexpr bigint<n> field_char();

    /**
     * Initializes euler, s, t, t_minus_1_over_2, nqr, and nqr_to_t. These constants
     * must be set before calling sqrt(), either with this function or manually.
     */
    static void init_tonelli_shanks_constants();

    /** If base field, is the identity. */
    T Frobenius_map(unsigned long power) const;

    /* Functions common to all finite fields */

#ifdef PROFILE_OP_COUNTS // NOTE: op counts are affected when you exponentiate with ^
    static long long add_cnt;
    static long long sub_cnt;
    static long long mul_cnt;
    static long long sqr_cnt;
    static long long inv_cnt;
#endif

    virtual T& operator+=(const T& other) = 0;
    virtual T& operator-=(const T& other) = 0;
    virtual T& operator*=(const T& other) = 0;
    virtual T& operator^=(const unsigned long pow) = 0;
    template<mp_size_t m>
    virtual T& operator^=(const bigint<m> &pow) = 0;

    virtual T& square() = 0; 
    virtual T& invert() = 0;

    virtual T operator+(const T& other) const;
    virtual T operator-(const T& other) const;
    virtual T operator*(const T& other) const;
    virtual T operator^(const unsigned long pow) const;
    template<mp_size_t m>
    virtual T operator^(const bigint<m> &pow) const;
    virtual T operator-() const = 0;

    virtual T squared() const;
    virtual T inverse() const;
    /** HAS TO BE A SQUARE (else does not terminate). */
    virtual T sqrt() const = 0;

    bool operator==(const T& other) const = 0;
    bool operator!=(const T& other) const = 0;
    bool is_zero() const = 0;

    void print() const = 0;

    void randomize() = 0;
    void clear() = 0;

    /* The static functions should be defined in field classes, but are static so they
       can't be inherited. */
    static T zero();
    static T one();
    static T random_element();
    /** Equals 1 for prime field Fp. */
    static constexpr std::size_t extension_degree();

    /**
     * Returns the constituent bits in 64 bit words, in little-endian order.
     * Only the right-most ceil_size_in_bits() bits are used; other bits are 0.
     */
    std::vector<uint64_t> to_words() const = 0;
    /**
     * Creates a field element from the given bits in 64 bit words, in little-endian order.
     * Only the right-most ceil_size_in_bits() bits are used; other bits are ignored.
     */
    static T from_words(std::vector<uint64_t> words);
    static std::size_t ceil_size_in_bits();
    static std::size_t floor_size_in_bits();

    // the following should be defined as well but can't be inherited
    friend std::ostream& operator<<(std::ostream &out, const T &p);
    friend std::istream& operator>>(std::istream &in, T &p);
};

} // libff
