/**@file
 *****************************************************************************
 Declaration of GF(2^32) finite field.
 *****************************************************************************
 * @author     This file is part of libff (see AUTHORS), migrated from libiop
 * @copyright  MIT license (see LICENSE file)
 *****************************************************************************/
#ifndef LIBFF_ALGEBRA_GF32_HPP_
#define LIBFF_ALGEBRA_GF32_HPP_

#include <cstddef>
#include <cstdint>
#include <libff/algebra/field_utils/bigint.hpp>
#include <vector>

namespace libff {

/* gf32 implements the field GF(2)/[x^32 + x^22 + x^2 + x^1 + 1].
   Elements are represented internally with a single uint32 */
class gf32 {
public:
#ifdef PROFILE_OP_COUNTS // NOTE: op counts are affected when you exponentiate with ^
    static long long add_cnt;
    static long long sub_cnt;
    static long long mul_cnt;
    static long long sqr_cnt;
    static long long inv_cnt;
#endif
    // x^32 + x^22 + x^2 + x^1 + 1
    static const constexpr uint64_t modulus_ = 0b10000000000000000000111;
    static const constexpr uint64_t num_bits = 32;

    explicit gf32();
    explicit gf32(const uint32_t value);

    gf32& operator+=(const gf32 &other);
    gf32& operator-=(const gf32 &other);
    gf32& operator*=(const gf32 &other);
    gf32& operator^=(const unsigned long pow);
    template<mp_size_t m>
    gf32& operator^=(const bigint<m> &pow);

    gf32& square();
    gf32& invert();

    gf32 operator+(const gf32 &other) const;
    gf32 operator-(const gf32 &other) const;
    gf32 operator-() const;
    gf32 operator*(const gf32 &other) const;
    gf32 operator^(const unsigned long pow) const;
    template<mp_size_t m>
    gf32 operator^(const bigint<m> &pow) const;

    gf32 squared() const;
    gf32 inverse() const;
    gf32 sqrt() const;

    void randomize();
    void clear();

    bool operator==(const gf32 &other) const;
    bool operator!=(const gf32 &other) const;

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

    static gf32 random_element();

    static gf32 zero();
    static gf32 one();
    static gf32 multiplicative_generator; // generator of gf32^*

    static std::size_t ceil_size_in_bits() { return num_bits; }
    static std::size_t floor_size_in_bits() { return num_bits; }
    static constexpr std::size_t extension_degree() { return 32; }
    template<mp_size_t n>
    static constexpr bigint<n> field_char() { return bigint<n>(2); }

    friend std::ostream& operator<<(std::ostream &out, const gf32 &p);
    friend std::istream& operator>>(std::istream &in, gf32 &p);
private:
    uint32_t value_;
};

#ifdef PROFILE_OP_COUNTS
long long gf32::add_cnt = 0;
long long gf32::sub_cnt = 0;
long long gf32::mul_cnt = 0;
long long gf32::sqr_cnt = 0;
long long gf32::inv_cnt = 0;
#endif

} // namespace libff
#include <libff/algebra/fields/binary/gf32.tcc>

#endif // #ifndef LIBFF_ALGEBRA_GF32_HPP_
