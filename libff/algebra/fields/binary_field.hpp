/** @file
 *****************************************************************************
 Declaration of common API for all finite fields in the binary/ directory.

 Currently NOT used by the fields in this library. This class is not actually
 the parent class of any field. All APIs are enforced through tests instead.

 The reason for this is to ensure high performance of all fields. This class
 exists as documentation for common API between fields.

 Includes fields F_{2^n} for some selected values of n. All of the binary
 entension fields must implement all functions declared in this class.
 *****************************************************************************
 * @author     This file is part of libff, developed by SCIPR Lab
 *             and contributors (see AUTHORS).
 * @copyright  MIT license (see LICENSE file)
 *****************************************************************************/
#include <libff/algebra/field_utils/bigint.hpp>
#include <vector>

namespace libff {

template<typename T>
class BinaryField;

/* The type parameter T is intended to be set to the child class
 * when this class is extended. For example,
 * class gf32 : public BinaryField<gf32> ...
 */
template<typename T>
class BinaryField {
public:
    /* Functions unique to binary fields */

    // TODO: add documentation about how moduli are represented.
    static const constexpr uint64_t modulus_;
    static const constexpr uint64_t num_bits;

    /** generator of gf2^n */
    static T multiplicative_generator;

    /** If extension field, returns the base field's characteristic. */
    template<mp_size_t n>
    static constexpr bigint<n> field_char() { return bigint<n>(2); }

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
    // Floor and ceil should return the same thing for binary fields.
    static std::size_t ceil_size_in_bits() { return num_bits; }
    static std::size_t floor_size_in_bits() { return num_bits; }

    // the following should be defined as well but can't be inherited
    friend std::ostream& operator<<(std::ostream &out, const T &p);
    friend std::istream& operator>>(std::istream &in, T &p);
};

} // libff
