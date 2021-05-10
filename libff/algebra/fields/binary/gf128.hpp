/**@file
 *****************************************************************************
 Declaration of GF(2^128) finite field.
 *****************************************************************************
 * @author     This file is part of libff (see AUTHORS), migrated from libiop
 * @copyright  MIT license (see LICENSE file)
 *****************************************************************************/

#ifndef LIBFF_ALGEBRA_GF128_HPP_
#define LIBFF_ALGEBRA_GF128_HPP_

#include <cstddef>
#include <cstdint>
#include <vector>
#include <libff/algebra/field_utils/bigint.hpp>

namespace libff {

/* gf128 implements the field GF(2)/(x^128 + x^7 + x^2 + x + 1).
   Elements are represented internally with two uint64s */
class gf128 {
public:
#ifdef PROFILE_OP_COUNTS // NOTE: op counts are affected when you exponentiate with ^
    static long long add_cnt;
    static long long sub_cnt;
    static long long mul_cnt;
    static long long sqr_cnt;
    static long long inv_cnt;
#endif
    // x^128 + x^7 + x^2 + x + 1
    static const constexpr uint64_t modulus_ = 0b10000111;
    static const constexpr uint64_t num_bits = 128;

    explicit gf128();
    /* we need a constructor that only initializes the low half of value_ to
       be able to do gf128(0) and gf128(1). */
    explicit gf128(const uint64_t value_low);
    explicit gf128(const uint64_t value_high, const uint64_t value_low);

    gf128& operator+=(const gf128 &other);
    gf128& operator-=(const gf128 &other);
    gf128& operator*=(const gf128 &other);
    gf128& operator^=(const unsigned long pow);
    template<mp_size_t m>
    gf128& operator^=(const bigint<m> &pow);
    
    gf128& square();
    gf128& invert();

    gf128 operator+(const gf128 &other) const;
    gf128 operator-(const gf128 &other) const;
    gf128 operator-() const;
    gf128 operator*(const gf128 &other) const;
    gf128 operator^(const unsigned long pow) const;
    template<mp_size_t m>
    gf128 operator^(const bigint<m> &pow) const;

    gf128 squared() const;
    gf128 inverse() const;
    gf128 sqrt() const;

    void randomize();
    void clear();

    bool operator==(const gf128 &other) const;
    bool operator!=(const gf128 &other) const;

    bool is_zero() const;

    void print() const;
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

    static gf128 random_element();

    static gf128 zero();
    static gf128 one();
    static gf128 multiplicative_generator; // generator of gf128^*

    static std::size_t ceil_size_in_bits() { return num_bits; }
    static std::size_t floor_size_in_bits() { return num_bits; }
    static constexpr std::size_t extension_degree() { return 128; }
    template<mp_size_t n>
    static constexpr bigint<n> field_char() { return bigint<n>(2); }

    friend std::ostream& operator<<(std::ostream &out, const gf128 &el);
    friend std::istream& operator>>(std::istream &in, gf128 &el);
private:
    /* little-endian */
    uint64_t value_[2];
};

#ifdef PROFILE_OP_COUNTS
long long gf128::add_cnt = 0;
long long gf128::sub_cnt = 0;
long long gf128::mul_cnt = 0;
long long gf128::sqr_cnt = 0;
long long gf128::inv_cnt = 0;
#endif

} // namespace libff
#include <libff/algebra/fields/binary/gf128.tcc>

#endif // namespace libff_ALGEBRA_GF128_HPP_
