/**@file
 *****************************************************************************
 Declaration of GF(2^256) finite field.
 *****************************************************************************
 * @author     This file is part of libff (see AUTHORS), migrated from libiop
 * @copyright  MIT license (see LICENSE file)
 *****************************************************************************/
#ifndef LIBFF_ALGEBRA_GF256_HPP_
#define LIBFF_ALGEBRA_GF256_HPP_

#include <cstddef>
#include <cstdint>
#include <vector>
#include <libff/algebra/field_utils/bigint.hpp>

namespace libff {

/* x^256 + x^10 + x^5 + x^2 + 1 */
/* gf256 implements the field GF(2)/(x^256 + x^10 + x^5 + x^2 + 1).
   Elements are represented internally with four uint64s */
class gf256 {
public:
#ifdef PROFILE_OP_COUNTS // NOTE: op counts are affected when you exponentiate with ^
    static long long add_cnt;
    static long long sub_cnt;
    static long long mul_cnt;
    static long long sqr_cnt;
    static long long inv_cnt;
#endif
    // x^256 + x^10 + x^5 + x^2 + 1
    static const constexpr uint64_t modulus_ = 0b10000100101;
    static const constexpr uint64_t num_bits = 256;

    explicit gf256();
    /* we need a constructor that only initializes the low 64 bits of value_ to
       be able to do gf256(0) and gf256(1). */
    explicit gf256(const uint64_t value_low);
    explicit gf256(const uint64_t value_high, const uint64_t value_midh,
                   const uint64_t value_midl, const uint64_t value_low);

    gf256& operator+=(const gf256 &other);
    gf256& operator-=(const gf256 &other);
    gf256& operator*=(const gf256 &other);
    gf256& operator^=(const unsigned long pow);
    template<mp_size_t m>
    gf256& operator^=(const bigint<m> &pow);

    gf256& square();
    gf256& invert();

    gf256 operator+(const gf256 &other) const;
    gf256 operator-(const gf256 &other) const;
    gf256 operator-() const;
    gf256 operator*(const gf256 &other) const;
    gf256 operator^(const unsigned long pow) const;
    template<mp_size_t m>
    gf256 operator^(const bigint<m> &pow) const;

    gf256 squared() const;
    gf256 inverse() const;
    gf256 sqrt() const;

    void randomize();
    void clear();

    bool operator==(const gf256 &other) const;
    bool operator!=(const gf256 &other) const;

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

    static gf256 random_element();

    static gf256 zero();
    static gf256 one();
    static gf256 multiplicative_generator; // generator of gf256^*

    static std::size_t ceil_size_in_bits() { return num_bits; }
    static std::size_t floor_size_in_bits() { return num_bits; }
    static constexpr std::size_t extension_degree() { return 256; }
    template<mp_size_t n>
    static constexpr bigint<n> field_char() { return bigint<n>(2); }

    friend std::ostream& operator<<(std::ostream &out, const gf256 &p);
    friend std::istream& operator>>(std::istream &in, gf256 &p);
private:
    /* little-endian */
    uint64_t value_[4];
};

#ifdef PROFILE_OP_COUNTS
long long gf256::add_cnt = 0;
long long gf256::sub_cnt = 0;
long long gf256::mul_cnt = 0;
long long gf256::sqr_cnt = 0;
long long gf256::inv_cnt = 0;
#endif

} // namespace libff
#include <libff/algebra/fields/binary/gf256.tcc>

#endif // namespace libff_ALGEBRA_GF256_HPP_
