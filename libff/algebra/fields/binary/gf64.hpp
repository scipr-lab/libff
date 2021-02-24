/**@file
 *****************************************************************************
 Declaration of GF(2^64) finite field.
 *****************************************************************************
 * @author     This file is part of libff (see AUTHORS), migrated from libiop
 * @copyright  MIT license (see LICENSE file)
 *****************************************************************************/
#ifndef LIBFF_ALGEBRA_GF64_HPP_
#define LIBFF_ALGEBRA_GF64_HPP_

#include <cstddef>
#include <cstdint>
#include <vector>
#include <libff/algebra/field_utils/bigint.hpp>

namespace libff {

/* gf64 implements the field GF(2)/[x^64 + x^4 + x^3 + x + 1].
   Elements are represented internally with a single uint64 */
class gf64 {
public:
#ifdef PROFILE_OP_COUNTS // NOTE: op counts are affected when you exponentiate with ^
    static long long add_cnt;
    static long long sub_cnt;
    static long long mul_cnt;
    static long long sqr_cnt;
    static long long inv_cnt;
#endif
    // x^64 + x^4 + x^3 + x + 1. The assembly code assumes that no term other
    // than x^64 is greater than x^31, to enable faster multiplication.
    static const constexpr uint64_t modulus_ = 0b11011;
    static const constexpr uint64_t num_bits = 64;

    explicit gf64();
    explicit gf64(const uint64_t value);

    gf64& operator+=(const gf64 &other);
    gf64& operator-=(const gf64 &other);
    gf64& operator*=(const gf64 &other);
    gf64& operator^=(const unsigned long pow);
    template<mp_size_t m>
    gf64& operator^=(const bigint<m> &pow);

    gf64& square();
    gf64& invert();

    gf64 operator+(const gf64 &other) const;
    gf64 operator-(const gf64 &other) const;
    gf64 operator-() const;
    gf64 operator*(const gf64 &other) const;
    gf64 operator^(const unsigned long pow) const;
    template<mp_size_t m>
    gf64 operator^(const bigint<m> &pow) const;

    gf64 squared() const;
    gf64 inverse() const;
    gf64 sqrt() const;

    void randomize();
    void clear();

    bool operator==(const gf64 &other) const;
    bool operator!=(const gf64 &other) const;

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

    static gf64 random_element();

    static gf64 zero();
    static gf64 one();
    static gf64 multiplicative_generator; // generator of gf64^*

    static std::size_t ceil_size_in_bits() { return num_bits; }
    static std::size_t floor_size_in_bits() { return num_bits; }
    static constexpr std::size_t extension_degree() { return 64; }
    template<mp_size_t n>
    static constexpr bigint<n> field_char() { return bigint<n>(2); }

    friend std::ostream& operator<<(std::ostream &out, const gf64 &p);
    friend std::istream& operator>>(std::istream &in, gf64 &p);
private:
    uint64_t value_;
};

#ifdef PROFILE_OP_COUNTS
long long gf64::add_cnt = 0;
long long gf64::sub_cnt = 0;
long long gf64::mul_cnt = 0;
long long gf64::sqr_cnt = 0;
long long gf64::inv_cnt = 0;
#endif

} // namespace libff
#include <libff/algebra/fields/binary/gf64.tcc>

#endif // namespace libff_ALGEBRA_GF64_HPP_
