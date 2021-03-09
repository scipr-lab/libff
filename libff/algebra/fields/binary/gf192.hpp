/**@file
 *****************************************************************************
 Declaration of GF(2^192) finite field.
 *****************************************************************************
 * @author     This file is part of libff (see AUTHORS), migrated from libiop
 * @copyright  MIT license (see LICENSE file)
 *****************************************************************************/

#ifndef LIBFF_ALGEBRA_GF192_HPP_
#define LIBFF_ALGEBRA_GF192_HPP_

#include <cstddef>
#include <cstdint>
#include <vector>
#include <libff/algebra/field_utils/bigint.hpp>

namespace libff {

/* gf192 implements the field GF(2)/(x^192 + x^7 + x^2 + x + 1).
   Elements are represented internally with three uint64s */
class gf192 {
public:
#ifdef PROFILE_OP_COUNTS // NOTE: op counts are affected when you exponentiate with ^
    static long long add_cnt;
    static long long sub_cnt;
    static long long mul_cnt;
    static long long sqr_cnt;
    static long long inv_cnt;
#endif
    // x^192 + x^7 + x^2 + x + 1
    static const constexpr uint64_t modulus_ = 0b10000111;
    static const constexpr uint64_t num_bits = 192;

    explicit gf192();
    /* we need a constructor that only initializes the low half of value_ to
       be able to do gf192(0) and gf192(1). */
    explicit gf192(const uint64_t value_low);
    explicit gf192(const uint64_t value_high, const uint64_t value_mid, const uint64_t value_low);

    gf192& operator+=(const gf192 &other);
    gf192& operator-=(const gf192 &other);
    gf192& operator*=(const gf192 &other);
    gf192& operator^=(const unsigned long pow);
    template<mp_size_t m>
    gf192& operator^=(const bigint<m> &pow);
    
    gf192& square();
    gf192& invert();

    gf192 operator+(const gf192 &other) const;
    gf192 operator-(const gf192 &other) const;
    gf192 operator-() const;
    gf192 operator*(const gf192 &other) const;
    gf192 operator^(const unsigned long pow) const;
    template<mp_size_t m>
    gf192 operator^(const bigint<m> &pow) const;

    gf192 squared() const;
    gf192 inverse() const;
    gf192 sqrt() const;
    /**
     * Returns the constituent bits in 64 bit words, in little-endian order.
     * Only the right-most ceil_size_in_bits() bits are used; other bits are 0.
     */
    std::vector<uint64_t> to_words() const;
    /**
     * Sets the field element from the given bits in 64 bit words, in little-endian order.
     * Only the right-most ceil_size_in_bits() bits are used; other bits are ignored.
     * Should always return true since the right-most bits are always valid.
     */
    bool from_words(std::vector<uint64_t> words);

    void randomize();
    void clear();

    bool operator==(const gf192 &other) const;
    bool operator!=(const gf192 &other) const;

    bool is_zero() const;

    void print() const;

    static gf192 random_element();

    static gf192 zero();
    static gf192 one();
    static gf192 multiplicative_generator; // generator of gf192^*

    static std::size_t ceil_size_in_bits() { return num_bits; }
    static std::size_t floor_size_in_bits() { return num_bits; }
    static constexpr std::size_t extension_degree() { return 192; }
    template<mp_size_t n>
    static constexpr bigint<n> field_char() { return bigint<n>(2); }

    friend std::ostream& operator<<(std::ostream &out, const gf192 &p);
    friend std::istream& operator>>(std::istream &in, gf192 &p);
private:
    /* little-endian */
    uint64_t value_[3];
};

#ifdef PROFILE_OP_COUNTS
long long gf192::add_cnt = 0;
long long gf192::sub_cnt = 0;
long long gf192::mul_cnt = 0;
long long gf192::sqr_cnt = 0;
long long gf192::inv_cnt = 0;
#endif

} // namespace libff
#include <libff/algebra/fields/binary/gf192.tcc>

#endif // namespace libff_ALGEBRA_GF192_HPP_
